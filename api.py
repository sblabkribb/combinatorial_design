"""
Combinatorial Design Microservice

FastAPI service that wraps combinatorial DNA assembly design logic.
Accepts Part sequences + Module design via JSON, runs GoldenGate/Gibson
assembly simulation, and returns results as JSON with base64-encoded GenBank.
"""

import base64
import io
import os
import pathlib
import re
import shutil
import tempfile
from itertools import product as iter_product
from typing import Literal, Optional

from fastapi import FastAPI, HTTPException
from pydantic import BaseModel

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Restriction import HincII, BsaI
from Bio.SeqUtils import MeltingTemp as mt

from pydna.all import pcr, Assembly, Dseqrecord, assembly_fragments
from pydna.all import read as pydna_read
from pydna.primer import Primer
from pydna.common_sub_strings import terminal_overlap
from primers import create

app = FastAPI(title="Combinatorial Design Service", version="1.0.0")

ASSETS_DIR = pathlib.Path(__file__).parent / "assets"
DEFAULT_VECTOR_PATH = ASSETS_DIR / "pUC19.gb"


def _gibson_terminal_overlap_limit() -> int:
    """Return the minimum terminal overlap length used by pydna assembly."""
    raw = os.getenv("GIBSON_TERMINAL_OVERLAP_LIMIT", "20").strip()
    try:
        value = int(raw)
    except ValueError:
        return 20
    return max(1, value)


def _rc(seq: str) -> str:
    """Return reverse complement of a DNA sequence."""
    comp = str.maketrans("ATGCNatgcn", "TACGNtacgn")
    return seq.translate(comp)[::-1]


def _extract_primer_tail(primer_seq: str, template_seq: str, min_anneal: int = 12) -> str:
    """Return the 5' tail of primer_seq that does not anneal to template_seq."""
    p = primer_seq.upper()
    t = template_seq.upper()
    search_region = t[:100 + len(p)]
    for tail_len in range(len(p) - min_anneal + 1):
        anneal = p[tail_len:]
        check = min(len(anneal), min_anneal + 5)
        if anneal[:check] in search_region:
            return p[:tail_len]
    return p  # all tail, no annealing found


def _build_gg_pcr_product(
    template: str,
    fwd_seq: str,
    rev_seq: str,
    min_anneal: int = 12,
) -> str:
    """
    Reconstruct the Golden Gate PCR product from primers + template.

    fwd primer: 5'-[BsaI][OHL][extra][annealing_to_template_5']-3'
    rev primer: 5'-[BsaI_rev][OHR][extra][annealing_to_template_3'_rc]-3'

    PCR product = fwd_full + template_middle + rc(rev_full)
    where template_middle = template[len(fwd_anneal):len(template)-len(rev_anneal)]

    The product contains 2 BsaI sites, so _find_bsai_insert (Case A) handles it
    correctly and preserves the extra nt between overhang and annealing.
    """
    t = template.upper()
    f = fwd_seq.upper()
    r = rev_seq.upper()

    fwd_tail = _extract_primer_tail(f, t, min_anneal)
    fwd_anneal_len = len(f) - len(fwd_tail)

    rev_tail = _extract_primer_tail(r, _rc(t), min_anneal)
    rev_anneal_len = len(r) - len(rev_tail)

    middle_start = fwd_anneal_len
    middle_end = len(t) - rev_anneal_len
    middle = t[middle_start:middle_end] if middle_end > middle_start else ""

    return f + middle + _rc(r)


# ── Pydantic Schemas ──────────────────────────────────────────────────────────

class PartInput(BaseModel):
    id: str
    kid: str = ""
    name: str
    type: str
    sequence: str
    direction: str = "forward"
    overhang_left: str
    overhang_right: str
    source_type: Literal["part", "device"] = "part"
    source_id: Optional[str] = None
    primer_ids: list[str] = []


class ModuleSlot(BaseModel):
    type: str = ""
    part_names: list[str]


class ModuleInput(BaseModel):
    id: str
    slots: list[ModuleSlot]


class PrimerResult(BaseModel):
    target: str
    direction: str
    sequence: str
    tm: float
    length: int


class PrimerInputData(BaseModel):
    id: str
    name: str
    sequence: str
    direction: str = "forward"


class PrimerPosition(BaseModel):
    primer_id: str
    name: str
    sequence: str
    direction: str


class FeatureResult(BaseModel):
    label: str
    start: int
    end: int
    strand: int


class CombinationResult(BaseModel):
    id: str
    sequence: str
    length: int
    genbank_b64: str
    primers: list[PrimerResult]
    features: list[FeatureResult]


class PhaseConfig(BaseModel):
    phase: int
    assembly_method: str
    module_indices: list[int]


class RunRequest(BaseModel):
    project_name: str
    assembly_method: str = "goldengate"
    vector_sequence: Optional[str] = None
    vector_mcs_feature: str = "MCS"
    vector_overhang_left: Optional[str] = None   # GoldenGate: 벡터 좌측 4nt 오버행
    vector_overhang_right: Optional[str] = None  # GoldenGate: 벡터 우측 4nt 오버행
    vector_overlap_left: Optional[str] = None    # Gibson: 벡터 좌측 접합 오버랩 (20+ nt)
    vector_overlap_right: Optional[str] = None   # Gibson: 벡터 우측 접합 오버랩 (20+ nt)
    parts: list[PartInput]
    modules: list[ModuleInput]
    primer_ids: Optional[list[str]] = None
    primer_inputs: Optional[list[PrimerInputData]] = None
    phases: Optional[list[PhaseConfig]] = None   # 3G Assembly용 (미래 확장)


class RunResponse(BaseModel):
    combinations: list[CombinationResult]
    total_combinations: int
    errors: list[str]
    primer_positions: list[PrimerPosition] = []


# ── Helper Functions ──────────────────────────────────────────────────────────

def get_positions_by_label(record, labels: list[str]) -> dict:
    """Find bounding box of all features matching any label in the list."""
    start_list, end_list = [], []
    for feature in record.features:
        for label in labels:
            if "label" in feature.qualifiers:
                feat_label = feature.qualifiers["label"]
                if isinstance(feat_label, list):
                    feat_label = feat_label[0]
                if feat_label == label:
                    start_list.append(int(feature.location.start))
                    end_list.append(int(feature.location.end))
    if start_list and end_list:
        return {"start": min(start_list), "end": max(end_list)}
    return {"start": None, "end": None}


def record_to_genbank_b64(record: SeqRecord) -> str:
    """Serialize a SeqRecord to base64-encoded GenBank string."""
    buf = io.StringIO()
    SeqIO.write(record, buf, "genbank")
    return base64.b64encode(buf.getvalue().encode()).decode()


def extract_features(record) -> list[FeatureResult]:
    """Extract features from a record into serializable form."""
    results = []
    for f in record.features:
        label = ""
        if "label" in f.qualifiers:
            lv = f.qualifiers["label"]
            label = lv[0] if isinstance(lv, list) else lv
        results.append(FeatureResult(
            label=label,
            start=int(f.location.start),
            end=int(f.location.end),
            strand=f.location.strand if f.location.strand is not None else 1,
        ))
    return results


def _feature_labels(record) -> list[str]:
    labels: list[str] = []
    for feature in getattr(record, "features", []) or []:
        if "label" not in feature.qualifiers:
            continue
        value = feature.qualifiers["label"]
        label = value[0] if isinstance(value, list) else value
        labels.append(str(label))
    return labels


def _pick_goldengate_fragment(
    cut_fragments: list,
    part_label: str,
    overhang_left: str,
    overhang_right: str,
):
    """Pick the digested fragment that best matches the intended insert."""
    required = {part_label, overhang_left.upper(), overhang_right.upper()}
    part_only = []
    exact = []

    for frag in cut_fragments:
        labels = set(_feature_labels(frag))
        if required.issubset(labels):
            exact.append(frag)
        elif part_label in labels:
            part_only.append(frag)

    if exact:
        return max(exact, key=lambda x: len(x))
    if part_only:
        return max(part_only, key=lambda x: len(x))
    return max(cut_fragments, key=lambda x: len(x))


def _format_goldengate_declared_chain(gb_files: list[str], vector_used: bool) -> str:
    chain: list[str] = []
    if gb_files:
        first = gb_files[0].replace("-withvector.gb", "").split("-")
        last = gb_files[-1].replace("-withvector.gb", "").split("-")
        if vector_used and len(first) >= 2 and len(last) >= 2:
            chain.append(f"vector_backbone({last[-1]}->{first[0]})")

    for gbfile in gb_files:
        stem = gbfile.replace("-withvector.gb", "")
        parts = stem.split("-")
        if len(parts) >= 3:
            name = "-".join(parts[1:-1])
            chain.append(f"{name}({parts[0]}->{parts[-1]})")
        else:
            chain.append(stem)

    return " -> ".join(chain)


def _format_goldengate_fragment_summary(fragments: list) -> str:
    summary: list[str] = []
    for idx, frag in enumerate(fragments):
        try:
            labels = _feature_labels(frag)
            core_labels = [
                label for label in labels
                if label != "BsaI" and not re.fullmatch(r"[ACGTacgt]{3,8}", label)
            ]
            label = core_labels[0] if core_labels else (labels[0] if labels else f"frag{idx}")
            ds = frag.seq if hasattr(frag, "seq") else frag
            try:
                l_type, l_seq = ds.five_prime_end()
                r_type, r_seq = ds.three_prime_end()
                ends = f"L={l_type}({l_seq or 'blunt'}) R={r_type}({r_seq or 'blunt'})"
            except Exception:
                ends = f"ovhg={getattr(ds, 'ovhg', '?')}"
            summary.append(f"{label}({ends}, len={len(frag)})")
        except Exception as exc:
            summary.append(f"frag{idx}(summary_error={exc})")
    return "\n  ".join(summary)


def _reannotate_goldengate_assembly(
    assembled,
    gb_files: list[str],
    withvector_path: str,
    parts: list,
) -> None:
    """
    Restore part-name features that may be lost during pydna assembly.

    Strategy (tried in order):
    1. Read withvector .gb file with BioPython, find the part-labeled feature,
       use that exact sequence for the search.
    2. Fallback: use orig_part.sequence (template) from the request.
    In both cases, also try reverse complement if forward search fails.
    """
    assembled_seq = str(assembled.seq).upper()
    seq_len = len(assembled_seq)
    part_lookup = {p.name.replace(" ", "_"): p for p in parts}
    existing = {f.qualifiers.get("label", [""])[0] for f in assembled.features}

    for gbfile in gb_files:
        stem = gbfile.replace("-withvector.gb", "")
        stem_parts = stem.split("-")
        if len(stem_parts) < 3:
            continue
        part_name = "-".join(stem_parts[1:-1])
        if part_name in existing:
            continue

        orig_part = part_lookup.get(part_name)

        # --- Collect candidate sequences to search for ---
        candidates: list[str] = []

        # 1. Sequence from withvector file feature (preserves primer extra nt)
        try:
            bio_rec = SeqIO.read(os.path.join(withvector_path, gbfile), "genbank")
            for f in bio_rec.features:
                lbl = f.qualifiers.get("label", [""])[0]
                if lbl == part_name:
                    s = str(bio_rec.seq[int(f.location.start):int(f.location.end)]).upper()
                    if len(s) >= 4:
                        candidates.append(s)
                    break
        except Exception:
            pass

        # 2. Fallback: template sequence from request
        if orig_part and orig_part.sequence:
            seq = orig_part.sequence.upper()
            if orig_part.direction == "reverse":
                seq = _rc(seq)
            candidates.append(seq)

        if not candidates:
            continue

        # --- Search assembled sequence ---
        found_pos: int = -1
        found_seq: str = ""
        found_strand: int = 1

        doubled = assembled_seq + assembled_seq  # for wrap-around

        for cand in candidates:
            for strand, search_seq in [(1, cand), (-1, _rc(cand))]:
                pos = assembled_seq.find(search_seq)
                if pos == -1:
                    # try wrap-around
                    pos = doubled.find(search_seq)
                    if pos != -1 and pos < seq_len:
                        pass  # valid wrap-around start
                    else:
                        pos = -1
                if pos != -1:
                    found_pos = pos
                    found_seq = search_seq
                    found_strand = strand
                    break
            if found_pos != -1:
                break

        if found_pos == -1:
            continue

        feat_type = (orig_part.type if orig_part and orig_part.type else "misc_feature")
        end_pos = found_pos + len(found_seq)
        if end_pos > seq_len:
            # wrap-around feature — use CompoundLocation
            from Bio.SeqFeature import CompoundLocation
            loc = CompoundLocation([
                FeatureLocation(found_pos, seq_len, found_strand),
                FeatureLocation(0, end_pos - seq_len, found_strand),
            ])
        else:
            loc = FeatureLocation(found_pos, end_pos, found_strand)

        assembled.features.append(
            SeqFeature(loc, type=feat_type, qualifiers={"label": [part_name]})
        )


def remove_duplicate_features(record):
    """Remove features that share identical start/end positions."""
    seen: set[tuple] = set()
    unique = []
    for f in record.features:
        key = (int(f.location.start), int(f.location.end))
        if key not in seen:
            seen.add(key)
            unique.append(f)
    record.features = unique
    return record


def fix_feature_location_posix(new_construct, fragment_list, insert_amp):
    """
    Correct feature location shift on Linux/POSIX.
    pydna does not properly update feature locations after PCR primer extension.
    """
    if os.name != "posix":
        return new_construct
    updated_primer = fragment_list[1].forward_primer.seq
    org_primer = insert_amp.forward_primer.seq
    shift_len = len(updated_primer) - len(org_primer)
    if shift_len == 0:
        return new_construct

    insert_pos = get_positions_by_label(new_construct, ["insert_forward"])
    if insert_pos["start"] is None:
        return new_construct

    for feature in new_construct.features:
        if feature.location.start >= insert_pos["start"]:
            strand = feature.location.strand
            feature.location = FeatureLocation(
                start=feature.location.start + shift_len,
                end=feature.location.end + shift_len,
            )
            feature.location.strand = strand
    return new_construct


# ── Assembly Pipeline ─────────────────────────────────────────────────────────

def create_insert_genbank_files(parts: list[PartInput], insert_path: str) -> None:
    """
    Generate insert GenBank files from Part data.
    Each file is named {overhang_left}-{part.name}-{overhang_right}-insert.gb
    """
    for part in parts:
        part_id = f"{part.overhang_left}-{part.name}-{part.overhang_right}"
        sequence = part.sequence.upper().strip()
        if part.direction == "reverse":
            sequence = str(Seq(sequence).reverse_complement())

        seq_record = SeqRecord(
            Seq(sequence),
            id=part_id[:16],
            name=part_id[:16],
            description=f"Synthetic construct of {part.name}",
        )
        seq_record.features = [
            SeqFeature(
                FeatureLocation(0, len(sequence), strand=1),
                type="misc_feature",
                qualifiers={"label": [part.name]},
            )
        ]
        seq_record.annotations = {"molecule_type": "DNA"}
        SeqIO.write(
            seq_record,
            os.path.join(insert_path, f"{part_id}-insert.gb"),
            "genbank",
        )


# ── GoldenGate fragment preparation ──────────────────────────────────────────

def _find_bsai_insert(
    seq: str, overhang_left: str, overhang_right: str
) -> tuple[str | None, str]:
    """
    Extract the BsaI-flanked region [GGTCTC][N][OHL]...[OHR][N][GAGACC]
    from a sequence that has exactly 2 BsaI sites (1 fwd + 1 rev).

    Returns (region, "") on success or (None, error_message) on failure.
    error_message is empty string when there are simply no/wrong BsaI sites
    (caller falls back to wrap_with_bsai silently) and non-empty for
    explicit problems like >2 sites or overhang mismatch.
    """
    s = seq.upper()
    oh_l = overhang_left.upper()
    oh_r = overhang_right.upper()

    fwd = [i for i in range(len(s) - 5) if s[i:i + 6] == "GGTCTC"]
    rev = [i for i in range(len(s) - 5) if s[i:i + 6] == "GAGACC"]
    total = len(fwd) + len(rev)

    if total > 2:
        return None, (
            f"BsaI site {total}개 발견 (최대 2개 허용). Assembly에서 제외됩니다."
        )

    if len(fwd) == 1 and len(rev) == 1:
        fwd_pos, rev_pos = fwd[0], rev[0]
        if fwd_pos < rev_pos:
            region = s[fwd_pos: rev_pos + 6]
        else:
            # Circular wrap: GAGACC appears before GGTCTC in linear representation
            region = s[fwd_pos:] + s[: rev_pos + 6]

        # Verify overhangs: GGTCTC + 1 spacer = 7 chars before OHL
        actual_ohl = region[7: 7 + len(oh_l)]
        # GAGACC = 6 chars at end, preceded by 1 spacer → OHR ends at len-7
        ohr_end = len(region) - 7
        actual_ohr = region[ohr_end - len(oh_r): ohr_end]

        if actual_ohl != oh_l:
            return None, (
                f"left overhang 불일치 (기대: {oh_l}, 실제: {actual_ohl}). "
                "overhang 파라미터 기준으로 BsaI site를 추가합니다."
            )
        if actual_ohr != oh_r:
            return None, (
                f"right overhang 불일치 (기대: {oh_r}, 실제: {actual_ohr}). "
                "overhang 파라미터 기준으로 BsaI site를 추가합니다."
            )
        return region, ""

    return None, ""  # 0, 1, or 2 same-direction sites → fall back to wrap


def _make_withvector_record(
    construct_seq: str,
    safe_name: str,
    overhang_left: str,
    overhang_right: str,
    ohl_start: int,
    ohr_end: int,
) -> Dseqrecord:
    """Create an annotated Dseqrecord for a BsaI-flanked construct."""
    from Bio.SeqFeature import SeqFeature, FeatureLocation as FL

    oh_l = overhang_left.upper()
    oh_r = overhang_right.upper()
    ohl_end = ohl_start + len(oh_l)
    ohr_start = ohr_end - len(oh_r)

    record = Dseqrecord(construct_seq)
    record.features = [
        SeqFeature(FL(0, 6, 1), type="misc_feature", qualifiers={"label": ["BsaI"]}),
        SeqFeature(FL(ohl_start, ohl_end, 1), type="misc_feature", qualifiers={"label": [oh_l]}),
        SeqFeature(FL(ohl_end, ohr_start, 1), type="misc_feature", qualifiers={"label": [safe_name]}),
        SeqFeature(FL(ohr_start, ohr_end, 1), type="misc_feature", qualifiers={"label": [oh_r]}),
        SeqFeature(FL(ohr_end + 1, len(construct_seq), -1), type="misc_feature", qualifiers={"label": ["BsaI"]}),
    ]
    record.annotations["molecule_type"] = "DNA"
    record.name = safe_name[:16]
    record.id = safe_name[:16]
    record.description = ""
    return record


def _wrap_with_bsai(
    seq: str, overhang_left: str, overhang_right: str, safe_name: str
) -> Dseqrecord:
    """
    Case B: wrap a bare insert with BsaI sites and overhang sequences.
    Construct: 5'-GGTCTCA[OHL][insert][OHR]TGAGACC-3'
    After BsaI digestion: [OHL][insert][OHR] with 4-nt 5' overhangs.
    """
    oh_l = overhang_left.upper()
    oh_r = overhang_right.upper()
    construct = "GGTCTCA" + oh_l + seq.upper() + oh_r + "TGAGACC"
    # OHL starts at position 7 (GGTCTCA = 7 chars)
    # OHR ends at len(construct) - 7 (TGAGACC = 7 chars)
    return _make_withvector_record(construct, safe_name, oh_l, oh_r, 7, len(construct) - 7)


def prepare_goldengate_fragments(
    parts: list,
    modules: list,
    withvector_path: str,
    primer_inputs_map: dict | None = None,
) -> list[str]:
    """
    Prepare -withvector.gb files for GoldenGate assembly without Gibson cloning.

    Case A — part has exactly 2 BsaI sites (1 fwd GGTCTC + 1 rev GAGACC):
      Extract [GGTCTC][N][OHL][insert][OHR][N][GAGACC] as a linear fragment.
      Verify overhangs match slot parameters.
      Overhang mismatch → fallback to Case B with a warning.

    Case B — 0 or 1 BsaI sites, or overhang mismatch:
      Wrap bare sequence: GGTCTCA[OHL][seq][OHR]TGAGACC

    Primer PCR mode — part has primer_ids:
      Reconstruct PCR product = fwd_primer + template_middle + rc(rev_primer).
      The product embeds BsaI sites from primer tails (→ falls into Case A),
      preserving extra nucleotides between the overhang and annealing region.

    Error → skip part and record error message:
      > 2 BsaI sites in a single part.
    """
    errors: list[str] = []
    primer_inputs_map = primer_inputs_map or {}
    # overhang_left/right lives on PartInput, not ModuleSlot
    part_lookup = {p.name: p for p in parts}

    for module in modules:
        for slot in module.slots:
            for part_name in slot.part_names:
                part = part_lookup.get(part_name)
                if not part:
                    errors.append(f"Part '{part_name}' not found in request")
                    continue

                template_seq = part.sequence.upper()

                # ── Primer PCR mode ──────────────────────────────────────────
                # If primers are associated with this part, reconstruct the
                # full PCR product so that BsaI-OHL-extra and BsaI-OHR-extra
                # from the primer tails are included in the assembled fragment.
                if part.primer_ids and primer_inputs_map:
                    part_primers = [
                        primer_inputs_map[pid]
                        for pid in part.primer_ids
                        if pid in primer_inputs_map
                    ]
                    fwd_primers = [p for p in part_primers if p.direction == "forward"]
                    rev_primers = [p for p in part_primers if p.direction == "reverse"]

                    if fwd_primers and rev_primers:
                        seq = _build_gg_pcr_product(
                            template_seq,
                            fwd_primers[0].sequence,
                            rev_primers[0].sequence,
                        )
                    elif fwd_primers:
                        # Only forward primer available — prepend tail only
                        fwd = fwd_primers[0].sequence.upper()
                        tail = _extract_primer_tail(fwd, template_seq)
                        seq = fwd + template_seq[len(fwd) - len(tail):]
                    elif rev_primers:
                        # Only reverse primer available — append rc tail only
                        rev = rev_primers[0].sequence.upper()
                        tail = _extract_primer_tail(rev, _rc(template_seq))
                        seq = template_seq[: len(template_seq) - (len(rev) - len(tail))] + _rc(rev)
                    else:
                        seq = template_seq
                else:
                    seq = template_seq
                safe_name = part.name.replace(" ", "_")
                oh_l = part.overhang_left
                oh_r = part.overhang_right
                fname = f"{oh_l}-{safe_name}-{oh_r}-withvector.gb"
                fpath = os.path.join(withvector_path, fname)

                n_bsai = seq.count("GGTCTC") + seq.count("GAGACC")

                if n_bsai > 2:
                    errors.append(
                        f"{part.name}: BsaI site {n_bsai}개 발견 (최대 2개 허용). "
                        "Assembly에서 제외됩니다."
                    )
                    continue

                record = None

                if n_bsai == 2:
                    region, err = _find_bsai_insert(seq, oh_l, oh_r)
                    if region:
                        record = _make_withvector_record(
                            region, safe_name, oh_l, oh_r, 7, len(region) - 7,
                        )
                    elif err:
                        errors.append(f"{part.name}: {err}")

                if record is None:
                    record = _wrap_with_bsai(seq, oh_l, oh_r, safe_name)

                    # 래핑 후 접합부에서 우발적 BsaI 사이트 생성 여부 검사
                    # 예: seq 3' 말단 + OHR 5' 가 결합해 GAGACC/GGTCTC 생성
                    full_construct = "GGTCTCA" + oh_l + seq + oh_r + "TGAGACC"
                    n_construct = full_construct.count("GGTCTC") + full_construct.count("GAGACC")
                    if n_construct != 2:
                        extra = n_construct - 2
                        # OHR 접합부 확인
                        ohr_junction = seq[-6:] + oh_r
                        ohl_junction = oh_l + seq[:6]
                        if "GAGACC" in ohr_junction or "GGTCTC" in ohr_junction:
                            junction_hint = f"OHR '{oh_r}' + 서열 3' 말단 '{seq[-6:]}'"
                        elif "GAGACC" in ohl_junction or "GGTCTC" in ohl_junction:
                            junction_hint = f"OHL '{oh_l}' + 서열 5' 말단 '{seq[:6]}'"
                        else:
                            junction_hint = f"OHL/OHR 접합부"
                        errors.append(
                            f"[경고] {part.name}: {junction_hint}에서 "
                            f"의도치 않은 BsaI 사이트 {extra}개 생성 — "
                            "BsaI 절단 시 잘못된 fragment 선택될 수 있음. 오버행 변경 권장."
                        )

                SeqIO.write(record, fpath, "genbank")

    return errors


def prepare_gibson_fragments(
    parts: list,
    modules: list,
    withvector_path: str,
) -> list[str]:
    """
    Prepare -withvector.gb files for Gibson assembly.

    Each fragment = overlap_left + part_sequence + overlap_right.
    overlap_left/right are stored in PartInput.overhang_left/overhang_right.

    Validates:
    - overlap length >= 20 nt (warns if shorter)
    - slot type vs part type mismatch (warns)

    File naming: {ol[:8]}-{part_name}-{or[:8]}-withvector.gb
    (8nt prefix used in filename; filter_combinations_by_overhang reuses same logic)
    """
    errors: list[str] = []
    min_overlap = _gibson_terminal_overlap_limit()
    part_lookup = {p.name: p for p in parts}

    for module in modules:
        for slot in module.slots:
            for part_name in slot.part_names:
                part = part_lookup.get(part_name)
                if not part:
                    errors.append(f"Part '{part_name}' not found in request")
                    continue

                ol = part.overhang_left
                or_ = part.overhang_right

                # Overlap length validation — 3-tier rules (Item 5)
                skip = False
                if len(ol) < 10:
                    errors.append(
                        f"[에러] {part.name}: overlap_left 길이 {len(ol)}nt — "
                        "10nt 미만은 Gibson 조립 불가. 이 파트는 제외됩니다."
                    )
                    skip = True
                elif len(ol) < 20:
                    errors.append(
                        f"[경고] {part.name}: overlap_left 길이 {len(ol)}nt — "
                        f"Gibson은 20nt 이상 권장 (현재 엔진 최소 {min_overlap}nt, 10–19nt는 조립 불안정)"
                    )
                if len(or_) < 10:
                    errors.append(
                        f"[에러] {part.name}: overlap_right 길이 {len(or_)}nt — "
                        "10nt 미만은 Gibson 조립 불가. 이 파트는 제외됩니다."
                    )
                    skip = True
                elif len(or_) < 20:
                    errors.append(
                        f"[경고] {part.name}: overlap_right 길이 {len(or_)}nt — "
                        f"Gibson은 20nt 이상 권장 (현재 엔진 최소 {min_overlap}nt, 10–19nt는 조립 불안정)"
                    )
                if skip:
                    continue

                # Full fragment: overlap_left + core + overlap_right
                core_seq = part.sequence.upper()
                full_seq = ol.upper() + core_seq + or_.upper()
                safe_name = part.name.replace(" ", "_")

                # Use first 8nt of each overlap as filename key
                # (filter_combinations_by_overhang uses split("-")[0] and split("-")[-1])
                ol_key = ol[:8].upper()
                or_key = or_[:8].upper()
                fname = f"{ol_key}-{safe_name}-{or_key}-withvector.gb"
                fpath = os.path.join(withvector_path, fname)

                record = SeqRecord(
                    Seq(full_seq),
                    id=safe_name,
                    name=safe_name[:16],
                    description=f"Gibson fragment: {safe_name}",
                    annotations={"molecule_type": "DNA"},
                )
                ol_len = len(ol)
                core_len = len(core_seq)
                record.features = [
                    SeqFeature(
                        FeatureLocation(0, ol_len),
                        type="misc_feature",
                        qualifiers={"label": [ol_key]},
                    ),
                    SeqFeature(
                        FeatureLocation(ol_len, ol_len + core_len),
                        type="misc_feature",
                        qualifiers={"label": [safe_name]},
                    ),
                    SeqFeature(
                        FeatureLocation(ol_len + core_len, len(full_seq)),
                        type="misc_feature",
                        qualifiers={"label": [or_key]},
                    ),
                ]
                SeqIO.write(record, fpath, "genbank")

    return errors


def assemble_gibson(
    filtered: dict[str, tuple],
    withvector_path: str,
    module_insert_path: str,
    vector_path: str | None = None,
    vector_overlap_left: str | None = None,
    vector_overlap_right: str | None = None,
) -> tuple[list[CombinationResult], list[str]]:
    """
    Run Gibson assembly for each filtered combination.

    Fragments already contain overlap regions on both ends
    (prepared by prepare_gibson_fragments).

    When vector_path is provided:
    - Wrap vector: overlap_right_last + vector_seq + overlap_left_first
    - Assemble circular (assemble_circular)
    When no vector:
    - Assemble linear (assemble_linear)

    Uses terminal_overlap with a configurable minimum overlap length.
    """
    results: list[CombinationResult] = []
    errors: list[str] = []
    min_overlap = _gibson_terminal_overlap_limit()

    # Read vector sequence once
    vector_seq: str | None = None
    if vector_path:
        try:
            vec_rec = pydna_read(vector_path)
            vector_seq = str(vec_rec.seq).upper().replace("N", "A")
        except Exception as exc:
            errors.append(f"Vector read error: {exc}")

    for combo_key, gb_files in filtered.items():
        fragments: list = []
        failed = False

        # Collect part fragments
        for gbfile in gb_files:
            try:
                record = SeqIO.read(os.path.join(withvector_path, gbfile), "genbank")
                frag = Dseqrecord(str(record.seq).upper(), linear=True)
                # Carry forward only the core part feature (not overlap labels)
                _OHG_RE = re.compile(r'^[ACGTacgt]{3,8}$')
                frag.features = [
                    f for f in record.features
                    if not _OHG_RE.match(f.qualifiers.get("label", [""])[0])
                ]
                fragments.append(frag)
            except Exception as exc:
                errors.append(f"{combo_key}/{gbfile}: {type(exc).__name__}: {exc}")
                failed = True
                break

        if failed or not fragments:
            continue

        try:
            do_circular = False
            if vector_seq:
                # Extract terminal overlaps from the combination filenames
                first_stem = gb_files[0].replace("-withvector.gb", "")
                last_stem  = gb_files[-1].replace("-withvector.gb", "")
                ol_first_key = first_stem.split("-")[0]   # 8nt prefix of first part's OL
                or_last_key  = last_stem.split("-")[-1]   # 8nt prefix of last part's OR

                # Validate vector_overlap_left/right if specified
                overhang_mismatch = False
                if vector_overlap_left and not vector_overlap_left.upper().startswith(ol_first_key):
                    errors.append(
                        f"{combo_key}: 벡터 좌측 오버랩 불일치 — "
                        f"벡터 요구 '{vector_overlap_left[:8]}', 첫 번째 파트 OL '{ol_first_key}'. "
                        "조합이 제외됩니다."
                    )
                    overhang_mismatch = True
                if vector_overlap_right and not vector_overlap_right.upper().startswith(or_last_key):
                    errors.append(
                        f"{combo_key}: 벡터 우측 오버랩 불일치 — "
                        f"벡터 요구 '{vector_overlap_right[:8]}', 마지막 파트 OR '{or_last_key}'. "
                        "조합이 제외됩니다."
                    )
                    overhang_mismatch = True

                if overhang_mismatch:
                    failed = True
                else:
                    # Vector fragment: OR_last + vector + OL_first
                    # (so terminal_overlap finds: fragments[-1] end matches vec start,
                    #  vec end matches fragments[0] start)
                    ol_seq = vector_overlap_left.upper() if vector_overlap_left else ""
                    or_seq = vector_overlap_right.upper() if vector_overlap_right else ""
                    vec_full = or_seq + vector_seq + ol_seq
                    vec_frag = Dseqrecord(vec_full, linear=True)
                    vec_frag.features = [
                        SeqFeature(
                            FeatureLocation(len(or_seq), len(or_seq) + len(vector_seq)),
                            type="misc_feature",
                            qualifiers={"label": ["vector_backbone"]},
                        )
                    ]
                    fragments = [vec_frag] + fragments
                    do_circular = True

            if failed:
                continue

            asm = Assembly(fragments, algorithm=terminal_overlap, limit=min_overlap)
            if do_circular:
                candidate = asm.assemble_circular()
                mode = "circular"
            else:
                candidate = asm.assemble_linear()
                mode = "linear"

            if not candidate:
                errors.append(f"{combo_key}: {mode} assembly produced no candidates")
                continue

            assembled = remove_duplicate_features(candidate[0])

            # Strip assembly-artifact features
            _OHG_RE2 = re.compile(r'^[ACGTacgt]{3,8}$')
            seq_len = len(assembled.seq)
            assembled.features = [
                f for f in assembled.features
                if not (
                    f.qualifiers.get("label", [""])[0] in ("vector_backbone",)
                    or _OHG_RE2.match(f.qualifiers.get("label", [""])[0])
                    or int(f.location.start) >= int(f.location.end)
                    or int(f.location.end) > seq_len
                )
            ]

            assembled.write(os.path.join(module_insert_path, f"{combo_key}.gb"))

            bio_record = SeqRecord(
                Seq(str(assembled.seq)),
                id=combo_key[:16],
                name=combo_key[:16],
                description=combo_key,
                features=assembled.features,
                annotations={"molecule_type": "DNA"},
            )

            results.append(CombinationResult(
                id=combo_key,
                sequence=str(assembled.seq),
                length=len(assembled.seq),
                genbank_b64=record_to_genbank_b64(bio_record),
                primers=[],   # Gibson assembly는 별도 primer 설계 불필요 (overlap이 primer 역할)
                features=extract_features(assembled),
            ))

        except Exception as exc:
            errors.append(f"{combo_key}: assembly error: {exc}")

    return results, errors


def clone_parts_into_vector(
    insert_path: str,
    withvector_path: str,
    vector_path: str,
    mcs_feature: str = "MCS",
) -> tuple[list[dict], list[str]]:
    """
    Clone each insert into a vector using Gibson assembly.
    Adapted from part_preparation.part_withvector_gibson().
    """
    import glob

    puc19 = pydna_read(vector_path)
    if not puc19.circular:
        puc19 = puc19.looped()

    mcs_pos = get_positions_by_label(puc19, [mcs_feature])
    if mcs_pos["start"] is None:
        half = len(puc19) // 2
        vector_frag1 = puc19[:half]
        vector_frag2 = puc19[half:]
    else:
        vector_frag1 = puc19[:mcs_pos["start"]]
        vector_frag2 = puc19[mcs_pos["end"]:]

    try:
        tmp_vector = puc19.linearize(HincII)
    except TypeError:
        # HincII doesn't cut this vector — use MCS/midpoint linearization as fallback
        tmp_vector = vector_frag2 + vector_frag1

    _, vrev1 = create(str(vector_frag1.seq).upper().replace('N', 'A'))
    vfwd2, _ = create(str(vector_frag2.seq).upper().replace('N', 'A'))

    vfwd = Primer(vrev1.seq)
    vrev = Primer(vfwd2.seq)
    vfwd.name = "vector_forward"
    vrev.name = "vector_reverse"
    vector_amp = pcr(vfwd, vrev, tmp_vector)

    primer_list: list[dict] = []
    errors: list[str] = []

    for fn in glob.glob(os.path.join(insert_path, "*-insert.gb")):
        try:
            insert = Dseqrecord(pydna_read(fn))
            fwd, rev = create(str(insert.seq))
            ifwd, irev = Primer(fwd.seq), Primer(rev.seq)

            # Guard against primers longer than insert
            while len(insert.seq) < len(ifwd) + len(irev):
                ifwd = Primer(ifwd.seq[:-1])
                irev = Primer(irev.seq[:-1])

            ifwd.name = "insert_forward"
            irev.name = "insert_reverse"

            insert_amp = pcr(ifwd, irev, insert)

            vector = Dseqrecord(vector_amp.seq)
            vector.features = vector_amp.features

            fragment_list = assembly_fragments((vector, insert_amp, vector))
            fragment_list = [vector_amp, fragment_list[1]]

            if not primer_list:
                primer_list += [
                    {"target": "vector", "direction": "forward",
                     "sequence": str(fragment_list[0].forward_primer.seq),
                     "tm": float(mt.Tm_NN(fragment_list[0].forward_primer.seq)),
                     "length": len(fragment_list[0].forward_primer.seq)},
                    {"target": "vector", "direction": "reverse",
                     "sequence": str(fragment_list[0].reverse_primer.seq),
                     "tm": float(mt.Tm_NN(fragment_list[0].reverse_primer.seq)),
                     "length": len(fragment_list[0].reverse_primer.seq)},
                ]

            asm = Assembly(fragment_list, algorithm=terminal_overlap)
            candidate = asm.assemble_circular()
            if not candidate:
                errors.append(f"Gibson assembly produced no candidates for {os.path.basename(fn)}")
                continue

            primer_list += [
                {"target": insert.id, "direction": "forward",
                 "sequence": str(fragment_list[1].forward_primer.seq),
                 "tm": float(mt.Tm_NN(fragment_list[1].forward_primer.seq)),
                 "length": len(fragment_list[1].forward_primer.seq)},
                {"target": insert.id, "direction": "reverse",
                 "sequence": str(fragment_list[1].reverse_primer.seq),
                 "tm": float(mt.Tm_NN(fragment_list[1].reverse_primer.seq)),
                 "length": len(fragment_list[1].reverse_primer.seq)},
            ]

            new_construct = candidate[0]
            new_construct = fix_feature_location_posix(new_construct, fragment_list, insert_amp)
            new_construct = remove_duplicate_features(new_construct)

            out_name = os.path.basename(fn).replace("-insert.gb", "-withvector.gb")
            new_construct.write(os.path.join(withvector_path, out_name))

        except Exception as exc:
            errors.append(f"Error cloning {os.path.basename(fn)}: {exc}")

    return primer_list, errors


def filter_combinations_by_overhang(
    modules: list[ModuleInput],
    withvector_files: list[str],
) -> dict[str, tuple]:
    """
    Generate all part combinations across ALL modules (sequential) and filter
    by overhang compatibility.

    All slots from all modules are treated as an ordered sequence. The Cartesian
    product is taken across the entire slot list, and the overhang chain is
    validated end-to-end (including inter-module junctions).

    Key scheme: "{all_module_ids_joined}_{part_names_joined}"

    Example:
        Module 1: [PromA|PromB] — [RBS1]
        Module 2: [HcaR|HbpR]
        → 4 combinations: PromA-RBS1-HcaR, PromA-RBS1-HbpR, ...
    """
    filtered: dict[str, tuple] = {}

    # Flatten all slots from all modules in order
    all_slots: list[ModuleSlot] = [slot for module in modules for slot in module.slots]
    module_prefix = "_".join(m.id for m in modules)

    if not all_slots:
        return filtered

    # For each slot, collect matching withvector files
    slots_files: list[list[str]] = []
    for slot in all_slots:
        matching = [wf for wf in withvector_files
                    for pname in slot.part_names if pname in wf]
        slots_files.append(matching)

    for combo in iter_product(*slots_files):
        if len(combo) < 2:
            if combo:
                fname = combo[0]
                key = fname.replace("-withvector.gb", "")
                filtered[f"{module_prefix}_{key}"] = combo
            continue

        # Validate full overhang chain (including inter-module junctions)
        overhang_right = [f.replace("-withvector.gb", "").split("-")[-1] for f in combo[:-1]]
        overhang_left  = [f.replace("-withvector.gb", "").split("-")[0]  for f in combo[1:]]

        if overhang_right == overhang_left:
            key = "-".join(
                "-".join(f.replace("-withvector.gb", "").split("-")[1:-1])
                for f in combo
            )
            filtered[f"{module_prefix}_{key}"] = combo

    return filtered


def assemble_goldengate(
    filtered: dict[str, tuple],
    withvector_path: str,
    module_insert_path: str,
    vector_path: str | None = None,
    vector_overhang_left: str | None = None,
    vector_overhang_right: str | None = None,
    parts: list | None = None,
) -> tuple[list[CombinationResult], list[str]]:
    """
    Run Golden Gate assembly for each filtered combination.
    Adapted from part_assembly.part_assembly_goldengate().

    When vector_path is provided, the vector backbone is wrapped with BsaI sites
    matching the terminal overhangs of the combination and included in a circular assembly.
    """
    results: list[CombinationResult] = []
    errors: list[str] = []
    parts = parts or []

    # Read vector sequence once (used for all combinations)
    vector_seq: str | None = None
    if vector_path:
        try:
            vec_rec = pydna_read(vector_path)
            vector_seq = str(vec_rec.seq).upper().replace("N", "A")
        except Exception as exc:
            errors.append(f"Vector read error: {exc}")

    for combo_key, gb_files in filtered.items():
        fragments = []
        primer_list: list[PrimerResult] = []
        failed = False

        for gbfile in gb_files:
            try:
                record = pydna_read(os.path.join(withvector_path, gbfile))

                # Labels to find: overhang codes + part name + "BsaI"
                labels = gbfile.replace("-withvector.gb", "").split("-") + ["BsaI"]
                pos = get_positions_by_label(record, labels)

                if pos["start"] is None:
                    errors.append(f"{combo_key}: feature labels not found in {gbfile}")
                    failed = True
                    break

                record_sub = record[pos["start"]:pos["end"]]
                record_sub.features = [
                    f for f in record.features
                    if (f.location.start >= pos["start"]
                        and f.location.end <= pos["end"]
                        and f.type != "primer_bind")
                ]

                # PCR-amplify the construct to obtain a clean dsDNA template.
                # Falls back to record_sub directly when the template is too short
                # (e.g. synthetic GoldenGate constructs generated by prepare_goldengate_fragments).
                ifwd, irev = None, None
                try:
                    fwd, rev = create(str(record_sub.seq).upper().replace('N', 'A'))
                    ifwd, irev = Primer(fwd.seq), Primer(rev.seq)
                    ifwd.name = "insert_forward"
                    irev.name = "insert_reverse"
                    insert_amp = pcr(ifwd, irev, record_sub.seq)
                    insert_amp.features = record_sub.features

                    # Normalize feature positions to zero
                    if insert_amp.features:
                        min_pos = min(
                            min(int(f.location.start), int(f.location.end))
                            for f in insert_amp.features
                        )
                        for f in insert_amp.features:
                            s = f.location.strand
                            f.location = FeatureLocation(
                                int(f.location.start) - min_pos,
                                int(f.location.end) - min_pos,
                                strand=s,
                            )
                except Exception:
                    # Template too short for PCR — use record_sub sequence directly as linear dsDNA
                    insert_amp = Dseqrecord(str(record_sub.seq), linear=True)
                    insert_amp.features = record_sub.features

                cut_fragments = insert_amp.cut(BsaI)
                # Need at least 2 fragments to confirm BsaI actually cut the construct.
                # cut() returns [self] (1 element) when no site is found — treat that as failure.
                if len(cut_fragments) < 2:
                    errors.append(
                        f"{combo_key}: no BsaI site (GGTCTC) found in {gbfile}. "
                        "The sequence may lack a BsaI recognition site."
                    )
                    failed = True
                    break

                stem = gbfile.replace("-withvector.gb", "")
                stem_parts = stem.split("-")
                expected_ohl = stem_parts[0] if len(stem_parts) >= 3 else ""
                expected_ohr = stem_parts[-1] if len(stem_parts) >= 3 else ""
                expected_name = "-".join(stem_parts[1:-1]) if len(stem_parts) >= 3 else stem

                selected_fragment = _pick_goldengate_fragment(
                    cut_fragments,
                    part_label=expected_name,
                    overhang_left=expected_ohl,
                    overhang_right=expected_ohr,
                )
                fragments.append(selected_fragment)
                if ifwd is not None and irev is not None:
                    primer_list += [
                        PrimerResult(
                            target=gbfile.replace("-withvector.gb", ""),
                            direction="forward",
                            sequence=str(ifwd.seq),
                            tm=float(mt.Tm_NN(ifwd.seq)),
                            length=len(ifwd.seq),
                        ),
                        PrimerResult(
                            target=gbfile.replace("-withvector.gb", ""),
                            direction="reverse",
                            sequence=str(irev.seq),
                            tm=float(mt.Tm_NN(irev.seq)),
                            length=len(irev.seq),
                        ),
                    ]

            except Exception as exc:
                errors.append(f"{combo_key}/{gbfile}: {type(exc).__name__}: {exc}")
                failed = True
                break

        if failed or not fragments:
            continue

        try:
            # Build vector backbone fragment for circular GoldenGate assembly.
            # Extract terminal overhangs from the ordered combination filenames:
            #   first file = OHL_first-name-OHR-withvector.gb  → OHL_first = split[0]
            #   last  file = OHL-name-OHR_last-withvector.gb   → OHR_last  = split[-1]
            do_circular = False
            if vector_seq:
                try:
                    first_stem = gb_files[0].replace("-withvector.gb", "")
                    last_stem  = gb_files[-1].replace("-withvector.gb", "")
                    ohl_first = first_stem.split("-")[0]
                    ohr_last  = last_stem.split("-")[-1]

                    # Q4: 벡터 오버행 요구사항 검증
                    # vector_overhang_left  = 벡터가 첫 번째 파트 왼쪽과 맞닿는 오버행 (= ohl_first와 일치해야 함)
                    # vector_overhang_right = 벡터가 마지막 파트 오른쪽과 맞닿는 오버행 (= ohr_last와 일치해야 함)
                    overhang_mismatch = False
                    if vector_overhang_left and ohl_first.upper() != vector_overhang_left.upper():
                        errors.append(
                            f"{combo_key}: 벡터 좌측 오버행 불일치 — "
                            f"벡터 요구 '{vector_overhang_left}', 첫 번째 파트 OHL '{ohl_first}'. "
                            "조합이 제외됩니다."
                        )
                        overhang_mismatch = True
                    if vector_overhang_right and ohr_last.upper() != vector_overhang_right.upper():
                        errors.append(
                            f"{combo_key}: 벡터 우측 오버행 불일치 — "
                            f"벡터 요구 '{vector_overhang_right}', 마지막 파트 OHR '{ohr_last}'. "
                            "조합이 제외됩니다."
                        )
                        overhang_mismatch = True
                    if overhang_mismatch:
                        failed = True
                    else:
                        # Silence internal BsaI sites in vector (GoldenGate requires
                        # BsaI-free backbone). Replace with 1-bp-variant that breaks
                        # recognition without changing backbone function in simulation.
                        n_internal = vector_seq.count("GGTCTC") + vector_seq.count("GAGACC")
                        if n_internal > 0:
                            vec_for_gg = (vector_seq
                                          .replace("GGTCTC", "GGTCTG")
                                          .replace("GAGACC", "GAGACG"))
                            errors.append(
                                f"{combo_key}: 벡터 내 BsaI 사이트 {n_internal}개 발견 — "
                                "시뮬레이션을 위해 자동 치환됨 (실제 실험 시 벡터 수정 필요)"
                            )
                        else:
                            vec_for_gg = vector_seq

                        # Wrap vector: GGTCTCA[OHR_last][vector][OHL_first]TGAGACC
                        # After BsaI cut the backbone has:
                        #   left  5' overhang = OHR_last  (ligates with last  part's right end)
                        #   right 5' overhang = OHL_first (ligates with first part's left  end)
                        vec_construct = "GGTCTCA" + ohr_last + vec_for_gg + ohl_first + "TGAGACC"
                        vec_dseq = Dseqrecord(vec_construct, linear=True)
                        vec_cut = vec_dseq.cut(BsaI)
                        if len(vec_cut) >= 2:
                            vec_frag = max(vec_cut, key=lambda x: len(x))
                            # Annotate the backbone so it's identifiable in the result
                            vec_frag.features = [
                                SeqFeature(
                                    FeatureLocation(0, len(vec_frag), 1),
                                    type="misc_feature",
                                    qualifiers={"label": ["vector_backbone"]},
                                )
                            ]
                            fragments = [vec_frag] + fragments
                            do_circular = True
                        else:
                            errors.append(f"{combo_key}: BsaI could not cut vector construct — assembling parts only")
                except Exception as vec_exc:
                    errors.append(f"{combo_key}: vector backbone prep failed: {vec_exc} — assembling parts only")

            if failed:
                continue

            asm = Assembly(fragments, algorithm=terminal_overlap, limit=4)
            if do_circular:
                candidate = asm.assemble_circular()
                mode = "circular"
            else:
                candidate = asm.assemble_linear()
                mode = "linear"
            if not candidate:
                errors.append(
                    f"{combo_key}: {mode} assembly produced no candidates.\n"
                    f"Declared chain: {_format_goldengate_declared_chain(gb_files, do_circular)}\n"
                    f"Actual fragments: {_format_goldengate_fragment_summary(fragments)}\n"
                    "Likely causes: internal BsaI site, incorrect cut-fragment selection, "
                    "or sticky-end orientation mismatch."
                )
                continue

            assembled = remove_duplicate_features(candidate[0])

            # Strip assembly-artifact features from the final plasmid:
            #  - BsaI recognition site markers
            #  - 4-nt overhang labels (e.g. "AAAA", "TCAC")
            #  - vector_backbone annotation added during assembly
            #  - any feature whose coordinates are out of range (circular wrap artifact)
            _OHG_RE = re.compile(r'^[ACGTacgt]{3,6}$')
            seq_len = len(assembled.seq)
            assembled.features = [
                f for f in assembled.features
                if not (
                    f.qualifiers.get("label", [""])[0] in ("BsaI", "vector_backbone")
                    or _OHG_RE.match(f.qualifiers.get("label", [""])[0])
                    or int(f.location.start) >= int(f.location.end)
                    or int(f.location.end) > seq_len
                )
            ]

            # Restore part features that may have been dropped by pydna
            # during restriction digestion + assembly
            _reannotate_goldengate_assembly(
                assembled, gb_files, withvector_path, parts
            )

            assembled.write(os.path.join(module_insert_path, f"{combo_key}.gb"))

            bio_record = SeqRecord(
                Seq(str(assembled.seq)),
                id=combo_key[:16],
                name=combo_key[:16],
                description=combo_key,
                features=assembled.features,
                annotations={"molecule_type": "DNA"},
            )

            results.append(CombinationResult(
                id=combo_key,
                sequence=str(assembled.seq),
                length=len(assembled.seq),
                genbank_b64=record_to_genbank_b64(bio_record),
                primers=primer_list,
                features=extract_features(assembled),
            ))

        except Exception as exc:
            errors.append(f"{combo_key}: assembly error: {exc}")

    return results, errors


# ── Endpoints ─────────────────────────────────────────────────────────────────

@app.get("/health")
async def health():
    return {"status": "ok"}


@app.post("/run", response_model=RunResponse)
async def run_assembly(request: RunRequest):
    """
    Run combinatorial assembly design.

    Steps:
    1. Generate insert GenBank files from Part sequences
    2. Clone inserts into vector (Gibson assembly)
    3. Filter part combinations by overhang compatibility
    4. Assemble valid combinations (GoldenGate / BsaI)
    5. Return assembled sequences, primers, and base64 GenBank files
    """
    tmp_dir = tempfile.mkdtemp(prefix="combdesign_")

    try:
        insert_path       = os.path.join(tmp_dir, "parts", "insert")
        withvector_path   = os.path.join(tmp_dir, "parts", "withvector")
        module_insert_path = os.path.join(tmp_dir, "modules", "inserts")
        for p in (insert_path, withvector_path, module_insert_path):
            os.makedirs(p, exist_ok=True)

        # Resolve vector
        if request.vector_sequence:
            vector_path = os.path.join(tmp_dir, "vector.gb")
            vec_record = SeqRecord(
                Seq(request.vector_sequence.upper()),
                id="custom_vector", name="custom_vector",
                description="Custom vector",
                annotations={"molecule_type": "DNA"},
            )
            SeqIO.write(vec_record, vector_path, "genbank")
        elif DEFAULT_VECTOR_PATH.exists():
            vector_path = str(DEFAULT_VECTOR_PATH)
        else:
            raise HTTPException(
                status_code=422,
                detail=(
                    "No vector provided and default pUC19 not found. "
                    "Provide vector_sequence or place pUC19.gb in assets/."
                ),
            )

        errors: list[str] = []

        # Step 1 & 2: prepare fragments per assembly method
        if request.assembly_method == "3g":
            return RunResponse(
                combinations=[],
                total_combinations=0,
                errors=["3G Assembly (GoldenGate + Gibson 결합)는 준비 중입니다."],
            )
        elif request.assembly_method == "goldengate":
            primer_inputs_map = {
                p.id: p for p in (request.primer_inputs or [])
            }
            errors.extend(prepare_goldengate_fragments(
                request.parts, request.modules, withvector_path,
                primer_inputs_map=primer_inputs_map,
            ))
        elif request.assembly_method == "gibson":
            errors.extend(prepare_gibson_fragments(
                request.parts, request.modules, withvector_path
            ))
        else:
            create_insert_genbank_files(request.parts, insert_path)
            _, clone_errors = clone_parts_into_vector(
                insert_path, withvector_path, vector_path,
                mcs_feature=request.vector_mcs_feature,
            )
            errors.extend(clone_errors)

        withvector_files = os.listdir(withvector_path)
        if not withvector_files:
            return RunResponse(
                combinations=[],
                total_combinations=0,
                errors=errors + [
                    "No withvector files generated. "
                    "Sequences may be too short for assembly or vector cloning failed."
                ],
            )

        # Step 3: filter combinations by overhang/overlap compatibility
        filtered = filter_combinations_by_overhang(request.modules, withvector_files)
        if not filtered:
            return RunResponse(
                combinations=[],
                total_combinations=0,
                errors=errors + [
                    "No valid combinations after overhang filtering. "
                    "Verify that adjacent parts share matching overhangs (GoldenGate: 4nt, Gibson: 20+nt prefix)."
                ],
            )

        # Step 4: assemble
        if request.assembly_method == "gibson":
            results, asm_errors = assemble_gibson(
                filtered, withvector_path, module_insert_path,
                vector_path=vector_path,
                vector_overlap_left=request.vector_overlap_left,
                vector_overlap_right=request.vector_overlap_right,
            )
        else:
            results, asm_errors = assemble_goldengate(
                filtered, withvector_path, module_insert_path,
                vector_path=vector_path,
                vector_overhang_left=request.vector_overhang_left,
                vector_overhang_right=request.vector_overhang_right,
                parts=request.parts,
            )
        errors.extend(asm_errors)

        # Build primer_positions from primer_inputs (coordinates computed on frontend)
        primer_positions: list[PrimerPosition] = []
        if request.primer_inputs:
            primer_input_map = {p.id: p for p in request.primer_inputs}
            for pid in (request.primer_ids or [p.id for p in request.primer_inputs]):
                if pid in primer_input_map:
                    p = primer_input_map[pid]
                    primer_positions.append(PrimerPosition(
                        primer_id=p.id, name=p.name,
                        sequence=p.sequence, direction=p.direction,
                    ))
                else:
                    errors.append(f"[경고] primer_id '{pid}' 가 primer_inputs에 없음 — 좌표 계산 불가")
        elif request.primer_ids:
            for pid in request.primer_ids:
                errors.append(f"[경고] primer_id '{pid}' 요청됨 — primer_inputs 미제공으로 좌표 계산 불가")

        return RunResponse(
            combinations=results,
            total_combinations=len(results),
            errors=errors,
            primer_positions=primer_positions,
        )

    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)
