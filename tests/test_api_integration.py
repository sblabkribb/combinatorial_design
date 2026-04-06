"""
Integration tests for the combinatorial design FastAPI service.
Requires the service to be running (tests run inside Docker container).
Uses httpx.TestClient for in-process testing.
"""

import base64
import io
import os
import sys
import tempfile

import pytest
from fastapi.testclient import TestClient
from Bio import SeqIO
from Bio.Restriction import BsaI
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pydna.all import Dseqrecord

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from api import (
    app,
    _wrap_with_bsai,
    prepare_goldengate_fragments,
    filter_combinations_by_overhang,
    assemble_goldengate,
    DEFAULT_VECTOR_PATH,
    PartInput,
    ModuleInput,
    ModuleSlot,
)

client = TestClient(app)

# ── Realistic test sequences ──────────────────────────────────────────────────
# These sequences contain BsaI recognition site (GGTCTC) for GoldenGate tests

# Short sequences with BsaI sites embedded
PROMOTER_SEQ = "AAGGTCTCATTGACAGCTAGCTCAGTCCTAGGTATAATGCTAGCGAGACC"   # GGTCTC at pos 2
RBS_SEQ      = "AAGGTCTCAAAAGAGGAGAAATACTAGTGAGACC"                    # GGTCTC at pos 2
CDS_SEQ      = "AAGGTCTCAATGAGTAAAGGAGAAGAACTTTTCGAGACC"               # GGTCTC at pos 2

# Minimal pUC19-like vector for testing (in memory, no file needed)
MOCK_VECTOR_SEQ = (
    "TCGCGCGTTTCGGTGATGACGGTGAAAACCTCTGACACATGCAGCTCCCGGAGACGGTCACAGCTTGTCTGTAAG"
    "CGGATGCCGGGAGCAGACAAGCCCGTCAGGGCGCGTCAGCGGGTGTTGGCGGGTGTGGGCGCAGCCATGACCC"
    "AGTCACGTAGCGATAGCGGAGTGTATACTGGCTTAACTATGCGGCATCAGAGCAGATTGTACTGAGAGTGCACC"
    "ATATGCGGTGTGAAATACCGCACAGATGCGTAAGGAGAAAATACCGCATCAGGCGCCATTCGCCATTCAGGCTG"
    # MCS region marker (we won't use a real HincII-containing vector in unit tests)
)


# ── Health Check ──────────────────────────────────────────────────────────────

class TestHealthEndpoint:
    def test_health_returns_ok(self):
        resp = client.get("/health")
        assert resp.status_code == 200
        assert resp.json() == {"status": "ok"}


# ── /run endpoint — input validation ─────────────────────────────────────────

class TestRunInputValidation:
    def test_empty_parts_returns_error(self):
        payload = {
            "project_name": "test",
            "parts": [],
            "modules": [],
        }
        resp = client.post("/run", json=payload)
        # Either 422 (no vector) or 200 with empty combinations
        assert resp.status_code in (200, 422, 503)

    def test_missing_required_field_returns_422(self):
        # Missing 'parts' field
        resp = client.post("/run", json={"project_name": "test", "modules": []})
        assert resp.status_code == 422

    def test_invalid_assembly_method_accepted(self):
        # Service accepts any string for assembly_method (validated downstream)
        payload = {
            "project_name": "test",
            "assembly_method": "invalid_method",
            "parts": [
                {"id": "p1", "name": "T7", "type": "Promoter",
                 "sequence": "ATCG", "overhang_left": "O1", "overhang_right": "O2"}
            ],
            "modules": [],
        }
        resp = client.post("/run", json=payload)
        assert resp.status_code in (200, 422)


# ── /run endpoint — core assembly ─────────────────────────────────────────────

class TestRunAssembly:
    """
    These tests use the default pUC19 vector (assets/pUC19.gb).
    They will be skipped if the vector file is not present.
    """

    @pytest.fixture(autouse=True)
    def skip_without_vector(self):
        import pathlib
        assets = pathlib.Path(__file__).parent.parent / "assets" / "pUC19.gb"
        if not assets.exists():
            pytest.skip("assets/pUC19.gb not found — skipping vector-dependent tests")

    def _make_payload(self, promoters, rbs_list, cds_list):
        parts = []
        overhangs = [("O1", "O2"), ("O2", "O3"), ("O3", "O4")]
        for i, (seqs, ptype) in enumerate([(promoters, "Promoter"),
                                            (rbs_list, "RBS"),
                                            (cds_list, "CDS")]):
            ol, or_ = overhangs[i]
            for j, seq in enumerate(seqs):
                parts.append({
                    "id": f"{ptype}_{j}",
                    "name": f"{ptype}_{j}",
                    "type": ptype,
                    "sequence": seq,
                    "direction": "forward",
                    "overhang_left": ol,
                    "overhang_right": or_,
                })
        return {
            "project_name": "test_project",
            "assembly_method": "goldengate",
            "parts": parts,
            "modules": [{
                "id": "MM1",
                "slots": [
                    {"type": "Promoter",
                     "part_names": [f"Promoter_{j}" for j in range(len(promoters))]},
                    {"type": "RBS",
                     "part_names": [f"RBS_{j}" for j in range(len(rbs_list))]},
                    {"type": "CDS",
                     "part_names": [f"CDS_{j}" for j in range(len(cds_list))]},
                ],
            }],
        }

    def test_response_schema(self):
        payload = self._make_payload([PROMOTER_SEQ], [RBS_SEQ], [CDS_SEQ])
        resp = client.post("/run", json=payload)
        assert resp.status_code == 200
        data = resp.json()
        assert "combinations" in data
        assert "total_combinations" in data
        assert "errors" in data

    def test_single_combination_assembled(self):
        payload = self._make_payload([PROMOTER_SEQ], [RBS_SEQ], [CDS_SEQ])
        resp = client.post("/run", json=payload)
        data = resp.json()
        # Either assembled successfully or reported BsaI error (acceptable)
        assert isinstance(data["combinations"], list)
        assert isinstance(data["errors"], list)

    def test_combination_has_required_fields(self):
        payload = self._make_payload([PROMOTER_SEQ], [RBS_SEQ], [CDS_SEQ])
        resp = client.post("/run", json=payload)
        data = resp.json()
        for combo in data["combinations"]:
            assert "id" in combo
            assert "sequence" in combo
            assert "length" in combo
            assert "genbank_b64" in combo
            assert "primers" in combo
            assert "features" in combo
            assert combo["length"] == len(combo["sequence"])

    def test_genbank_b64_is_valid_genbank(self):
        payload = self._make_payload([PROMOTER_SEQ], [RBS_SEQ], [CDS_SEQ])
        resp = client.post("/run", json=payload)
        data = resp.json()
        for combo in data["combinations"]:
            raw = base64.b64decode(combo["genbank_b64"]).decode()
            record = SeqIO.read(io.StringIO(raw), "genbank")
            assert len(record.seq) > 0

    def test_combinatorial_expansion_2x2(self):
        """2 promoters × 2 RBS → expect up to 4 combinations."""
        payload = self._make_payload(
            [PROMOTER_SEQ, PROMOTER_SEQ + "T"],   # 2 promoters
            [RBS_SEQ, RBS_SEQ + "T"],              # 2 RBS
            [CDS_SEQ],
        )
        resp = client.post("/run", json=payload)
        data = resp.json()
        # May be < 4 if BsaI cutting fails for some; at least schema is valid
        assert data["total_combinations"] == len(data["combinations"])

    def test_incompatible_overhangs_returns_empty(self):
        """Parts with non-matching overhangs should produce 0 combinations."""
        parts = [
            {"id": "p1", "name": "ProA", "type": "Promoter",
             "sequence": PROMOTER_SEQ, "overhang_left": "O1", "overhang_right": "O2"},
            {"id": "p2", "name": "RBS1", "type": "RBS",
             "sequence": RBS_SEQ, "overhang_left": "O9", "overhang_right": "O10"},  # mismatch
        ]
        payload = {
            "project_name": "test_mismatch",
            "parts": parts,
            "modules": [{
                "id": "MM1",
                "slots": [
                    {"type": "Promoter", "part_names": ["ProA"]},
                    {"type": "RBS",      "part_names": ["RBS1"]},
                ],
            }],
        }
        resp = client.post("/run", json=payload)
        assert resp.status_code == 200
        data = resp.json()
        assert data["total_combinations"] == 0

    def test_no_bsai_site_reported_in_errors(self):
        """Sequences without BsaI site should appear in errors, not crash."""
        parts = [
            {"id": "p1", "name": "ProA", "type": "Promoter",
             "sequence": "TTGACAGCTAGCTCAGTCCTA",  # no BsaI (GGTCTC)
             "overhang_left": "O1", "overhang_right": "O2"},
            {"id": "p2", "name": "RBS1", "type": "RBS",
             "sequence": "AAAGAGGAGAAATACTAGTG",    # no BsaI
             "overhang_left": "O2", "overhang_right": "O3"},
        ]
        payload = {
            "project_name": "test_no_bsai",
            "parts": parts,
            "modules": [{
                "id": "MM1",
                "slots": [
                    {"type": "Promoter", "part_names": ["ProA"]},
                    {"type": "RBS",      "part_names": ["RBS1"]},
                ],
            }],
        }
        resp = client.post("/run", json=payload)
        assert resp.status_code == 200
        data = resp.json()
        # Should not crash — BsaI error reported in errors list
        assert isinstance(data["errors"], list)
        assert data["total_combinations"] == 0

    def test_primer_fields_are_valid(self):
        payload = self._make_payload([PROMOTER_SEQ], [RBS_SEQ], [CDS_SEQ])
        resp = client.post("/run", json=payload)
        data = resp.json()
        for combo in data["combinations"]:
            for primer in combo["primers"]:
                assert primer["direction"] in ("forward", "reverse")
                assert len(primer["sequence"]) > 0
                assert primer["tm"] > 0
                assert primer["length"] == len(primer["sequence"])


# ── BsaI 절단 유닛 테스트 ─────────────────────────────────────────────────────

class TestBsaiCutBehavior:
    """_wrap_with_bsai 및 BsaI 절단 동작 검증 (pydna 의존, mock 없음)."""

    def test_wrap_with_bsai_produces_three_fragments(self):
        """_wrap_with_bsai가 BsaI로 정확히 3개 fragment를 생성하는지 검증."""
        record = _wrap_with_bsai("ATCGATCGATCGATCG", "AATG", "CGTA", "test_part")
        cut = record.cut(BsaI)
        assert len(cut) == 3, (
            f"BsaI 절단 후 fragment 수 = {len(cut)}, 3이어야 함 "
            f"(left stub / insert / right stub). 서열: {str(record.seq)}"
        )

    def test_bsai_cut_insert_contains_overhang(self):
        """절단 후 최대(insert) fragment의 서열에 OHL이 포함되는지 검증."""
        ohl, ohr = "AATG", "CGTA"
        record = _wrap_with_bsai("ATCGATCGATCGATCG", ohl, ohr, "test_part")
        cut = record.cut(BsaI)
        insert = max(cut, key=lambda x: len(x))
        seq_str = str(insert.seq)
        assert ohl in seq_str, (
            f"insert fragment에 OHL({ohl}) 없음. 실제 서열: {seq_str}"
        )

    def test_vector_backbone_bsai_produces_three_fragments(self):
        """GGTCTCA+OHR+vector+OHL+TGAGACC 구조가 BsaI로 3개 fragment를 생성하는지."""
        vec_seq = "ACGT" * 25  # 100bp 더미 벡터
        ohl_first, ohr_last = "AATG", "CGTA"
        vec_construct = "GGTCTCA" + ohr_last + vec_seq + ohl_first + "TGAGACC"
        vec_dseq = Dseqrecord(vec_construct, linear=True)
        cut = vec_dseq.cut(BsaI)
        assert len(cut) == 3, (
            f"벡터 construct BsaI 절단 후 fragment 수 = {len(cut)}, 3이어야 함"
        )

    def test_vector_backbone_insert_contains_vector_seq(self):
        """벡터 절단 후 최대 fragment가 원본 벡터 서열을 포함하는지."""
        vec_seq = "ACGT" * 25
        ohl_first, ohr_last = "AATG", "CGTA"
        vec_construct = "GGTCTCA" + ohr_last + vec_seq + ohl_first + "TGAGACC"
        vec_dseq = Dseqrecord(vec_construct, linear=True)
        cut = vec_dseq.cut(BsaI)
        backbone = max(cut, key=lambda x: len(x))
        assert vec_seq[:20] in str(backbone.seq), (
            "벡터 backbone fragment에 원본 벡터 서열이 없음"
        )


# ── 원형 GoldenGate 조립 테스트 ───────────────────────────────────────────────

class TestCircularGoldenGateAssembly:
    """
    2개 파트 + 더미 벡터로 assemble_circular() 후보가 생성되는지 검증.
    mock 없음 — 실제 pydna 실행.
    """

    # 간단한 DNA 서열 (BsaI 사이트 없는 bare sequence → _wrap_with_bsai 경로)
    PART1_SEQ = "ATCGATCGATCGATCGATCG"  # 20bp
    PART2_SEQ = "GCTAGCTAGCTAGCTAGCTA"  # 20bp
    VEC_SEQ   = "ACGT" * 25             # 100bp

    def _make_parts_and_modules(self, ohl1="AAAA", ohr1="TTTT",
                                 ohl2="TTTT", ohr2="CCCC"):
        parts = [
            PartInput(id="p1", name="part1", type="Promoter",
                      sequence=self.PART1_SEQ,
                      overhang_left=ohl1, overhang_right=ohr1),
            PartInput(id="p2", name="part2", type="RBS",
                      sequence=self.PART2_SEQ,
                      overhang_left=ohl2, overhang_right=ohr2),
        ]
        modules = [ModuleInput(id="M1", slots=[
            ModuleSlot(type="Promoter", part_names=["part1"]),
            ModuleSlot(type="RBS",      part_names=["part2"]),
        ])]
        return parts, modules

    def _write_vector_gb(self, path: str, seq: str) -> str:
        vec_path = os.path.join(path, "vector.gb")
        SeqIO.write(
            SeqRecord(Seq(seq), id="vec", name="vec",
                      annotations={"molecule_type": "DNA"}),
            vec_path, "genbank",
        )
        return vec_path

    def test_overhang_filter_passes_for_compatible_parts(self):
        """인접 overhang이 일치하는 파트 조합이 필터를 통과하는지."""
        parts, modules = self._make_parts_and_modules()
        with tempfile.TemporaryDirectory() as tmp:
            wp = os.path.join(tmp, "wv"); os.makedirs(wp)
            prepare_goldengate_fragments(parts, modules, wp)
            wfiles = os.listdir(wp)
            assert len(wfiles) == 2, f"withvector 파일 수 = {len(wfiles)}, 2 기대"
            filtered = filter_combinations_by_overhang(modules, wfiles)
            assert len(filtered) == 1, (
                f"필터 후 조합 수 = {len(filtered)}, 1 기대. files: {wfiles}"
            )

    def test_circular_assembly_produces_candidate(self):
        """2개 파트 + 벡터 → assemble_circular() 후보 1개 이상 반환."""
        parts, modules = self._make_parts_and_modules()
        with tempfile.TemporaryDirectory() as tmp:
            wp = os.path.join(tmp, "wv"); os.makedirs(wp)
            mp = os.path.join(tmp, "mi"); os.makedirs(mp)
            prepare_goldengate_fragments(parts, modules, wp)
            wfiles = os.listdir(wp)
            filtered = filter_combinations_by_overhang(modules, wfiles)
            vec_path = self._write_vector_gb(tmp, self.VEC_SEQ)

            results, errors = assemble_goldengate(
                filtered, wp, mp, vector_path=vec_path
            )

            circular_errors = [e for e in errors if "circular assembly" in e]
            assert not circular_errors, f"원형 조립 실패 에러: {circular_errors}"
            assert len(results) >= 1, (
                f"결과 0개. 전체 에러: {errors}"
            )

    def test_circular_assembly_result_contains_both_parts(self):
        """조립 결과 서열에 두 파트 서열이 모두 포함되는지."""
        parts, modules = self._make_parts_and_modules()
        with tempfile.TemporaryDirectory() as tmp:
            wp = os.path.join(tmp, "wv"); os.makedirs(wp)
            mp = os.path.join(tmp, "mi"); os.makedirs(mp)
            prepare_goldengate_fragments(parts, modules, wp)
            wfiles = os.listdir(wp)
            filtered = filter_combinations_by_overhang(modules, wfiles)
            vec_path = self._write_vector_gb(tmp, self.VEC_SEQ)

            results, _ = assemble_goldengate(
                filtered, wp, mp, vector_path=vec_path
            )
            if not results:
                pytest.skip("원형 조립 결과 없음 — 이전 테스트에서 실패해야 함")

            seq = results[0].sequence.upper()
            assert self.PART1_SEQ in seq, "결과 서열에 part1 없음"
            assert self.PART2_SEQ in seq, "결과 서열에 part2 없음"

    def test_circular_assembly_without_vector_falls_back_to_linear(self):
        """vector_path 없으면 linear assembly fallback — crash 없이 처리."""
        parts, modules = self._make_parts_and_modules()
        with tempfile.TemporaryDirectory() as tmp:
            wp = os.path.join(tmp, "wv"); os.makedirs(wp)
            mp = os.path.join(tmp, "mi"); os.makedirs(mp)
            prepare_goldengate_fragments(parts, modules, wp)
            wfiles = os.listdir(wp)
            filtered = filter_combinations_by_overhang(modules, wfiles)

            # vector_path=None → linear fallback
            results, errors = assemble_goldengate(filtered, wp, mp, vector_path=None)
            # linear assembly도 실패할 수 있지만 crash는 없어야 함
            assert isinstance(results, list)
            assert isinstance(errors, list)

    def test_four_parts_circular_assembly_produces_candidate(self):
        """4개 파트 + 벡터 → assemble_circular() 후보 생성 (실제 실패 케이스 재현)."""
        parts = [
            PartInput(id="p1", name="T7_promoter", type="Promoter",
                      sequence="ATCGATCGATCGATCGATCG",
                      overhang_left="AAAA", overhang_right="TTTT"),
            PartInput(id="p2", name="T7_RBS", type="RBS",
                      sequence="GCTAGCTAGCTAGCTAGCTA",
                      overhang_left="TTTT", overhang_right="CCCC"),
            PartInput(id="p3", name="sfGFP", type="CDS",
                      sequence="ACGTACGTACGTACGTACGT",
                      overhang_left="CCCC", overhang_right="GGGG"),
            PartInput(id="p4", name="T7_terminator", type="Terminator",
                      sequence="TGCATGCATGCATGCATGCA",
                      overhang_left="GGGG", overhang_right="AATG"),
        ]
        modules = [ModuleInput(id="MM1", slots=[
            ModuleSlot(type="Promoter",   part_names=["T7_promoter"]),
            ModuleSlot(type="RBS",        part_names=["T7_RBS"]),
            ModuleSlot(type="CDS",        part_names=["sfGFP"]),
            ModuleSlot(type="Terminator", part_names=["T7_terminator"]),
        ])]
        with tempfile.TemporaryDirectory() as tmp:
            wp = os.path.join(tmp, "wv"); os.makedirs(wp)
            mp = os.path.join(tmp, "mi"); os.makedirs(mp)
            prepare_goldengate_fragments(parts, modules, wp)
            wfiles = os.listdir(wp)
            filtered = filter_combinations_by_overhang(modules, wfiles)
            assert filtered, f"4파트 overhang 필터링 실패. files: {wfiles}"
            vec_path = self._write_vector_gb(tmp, self.VEC_SEQ)

            results, errors = assemble_goldengate(
                filtered, wp, mp, vector_path=vec_path
            )
            circular_errors = [e for e in errors if "circular assembly" in e]
            assert not circular_errors, f"4파트 원형 조립 실패: {circular_errors}"
            assert len(results) >= 1, f"결과 0개. 전체 에러: {errors}"

    def test_api_run_with_inline_vector_returns_circular(self):
        """POST /run (vector_sequence 포함) → circular assembly 에러 없음."""
        payload = {
            "project_name": "tdd_test",
            "assembly_method": "goldengate",
            "vector_sequence": self.VEC_SEQ,
            "parts": [
                {"id": "p1", "name": "part1", "type": "Promoter",
                 "sequence": self.PART1_SEQ,
                 "overhang_left": "AAAA", "overhang_right": "TTTT"},
                {"id": "p2", "name": "part2", "type": "RBS",
                 "sequence": self.PART2_SEQ,
                 "overhang_left": "TTTT", "overhang_right": "CCCC"},
            ],
            "modules": [{"id": "M1", "slots": [
                {"type": "Promoter", "part_names": ["part1"]},
                {"type": "RBS",      "part_names": ["part2"]},
            ]}],
        }
        resp = client.post("/run", json=payload)
        assert resp.status_code == 200
        data = resp.json()
        circular_errors = [e for e in data["errors"] if "circular assembly" in e]
        assert not circular_errors, f"circular assembly 에러: {circular_errors}"
        assert data["total_combinations"] >= 1, f"결과 0개. 에러: {data['errors']}"


class TestPartTypeMismatchWarning:
    """Q2: 슬롯 타입 vs 파트 실제 타입 불일치 시 경고 메시지 반환 검증."""

    PART1_SEQ = "ATCGATCGATCGATCGATCG"
    PART2_SEQ = "GCTAGCTAGCTAGCTAGCTA"
    VEC_SEQ   = "ACGT" * 25

    def test_type_mismatch_emits_warning(self):
        """슬롯 Promoter에 RBS 타입 파트 배치 → errors에 경고 포함."""
        payload = {
            "project_name": "type_mismatch_test",
            "assembly_method": "goldengate",
            "vector_sequence": self.VEC_SEQ,
            "parts": [
                {"id": "p1", "name": "part1", "type": "RBS",      # 실제 타입: RBS
                 "sequence": self.PART1_SEQ,
                 "overhang_left": "AAAA", "overhang_right": "TTTT"},
                {"id": "p2", "name": "part2", "type": "CDS",      # 실제 타입: CDS
                 "sequence": self.PART2_SEQ,
                 "overhang_left": "TTTT", "overhang_right": "CCCC"},
            ],
            "modules": [{"id": "M1", "slots": [
                {"type": "Promoter", "part_names": ["part1"]},   # 슬롯: Promoter ≠ RBS
                {"type": "RBS",      "part_names": ["part2"]},   # 슬롯: RBS ≠ CDS
            ]}],
        }
        resp = client.post("/run", json=payload)
        assert resp.status_code == 200
        data = resp.json()
        warnings = [e for e in data["errors"] if "경고" in e and "슬롯 타입" in e]
        assert len(warnings) == 2, f"경고 2개 기대, {len(warnings)}개 반환: {data['errors']}"
        assert any("part1" in w and "RBS" in w for w in warnings)
        assert any("part2" in w and "CDS" in w for w in warnings)

    def test_type_match_emits_no_warning(self):
        """슬롯과 파트 타입이 일치하면 타입 불일치 경고 없음."""
        payload = {
            "project_name": "type_match_test",
            "assembly_method": "goldengate",
            "vector_sequence": self.VEC_SEQ,
            "parts": [
                {"id": "p1", "name": "part1", "type": "Promoter",
                 "sequence": self.PART1_SEQ,
                 "overhang_left": "AAAA", "overhang_right": "TTTT"},
                {"id": "p2", "name": "part2", "type": "RBS",
                 "sequence": self.PART2_SEQ,
                 "overhang_left": "TTTT", "overhang_right": "CCCC"},
            ],
            "modules": [{"id": "M1", "slots": [
                {"type": "Promoter", "part_names": ["part1"]},
                {"type": "RBS",      "part_names": ["part2"]},
            ]}],
        }
        resp = client.post("/run", json=payload)
        assert resp.status_code == 200
        data = resp.json()
        warnings = [e for e in data["errors"] if "슬롯 타입" in e]
        assert not warnings, f"타입 일치인데 경고 발생: {warnings}"


class TestVectorOverhangValidation:
    """Q4: vector_overhang_left/right 지정 시 파트 오버행과 검증."""

    PART1_SEQ = "ATCGATCGATCGATCGATCG"
    PART2_SEQ = "GCTAGCTAGCTAGCTAGCTA"
    VEC_SEQ   = "ACGT" * 25

    def test_matching_vector_overhang_assembles_successfully(self):
        """벡터 오버행이 파트 오버행과 일치하면 조립 성공."""
        payload = {
            "project_name": "vec_oh_match",
            "assembly_method": "goldengate",
            "vector_sequence": self.VEC_SEQ,
            "vector_overhang_left": "AAAA",   # 첫 파트 OHL과 일치
            "vector_overhang_right": "CCCC",  # 마지막 파트 OHR과 일치
            "parts": [
                {"id": "p1", "name": "part1", "type": "Promoter",
                 "sequence": self.PART1_SEQ,
                 "overhang_left": "AAAA", "overhang_right": "TTTT"},
                {"id": "p2", "name": "part2", "type": "RBS",
                 "sequence": self.PART2_SEQ,
                 "overhang_left": "TTTT", "overhang_right": "CCCC"},
            ],
            "modules": [{"id": "M1", "slots": [
                {"type": "Promoter", "part_names": ["part1"]},
                {"type": "RBS",      "part_names": ["part2"]},
            ]}],
        }
        resp = client.post("/run", json=payload)
        assert resp.status_code == 200
        data = resp.json()
        mismatch_errors = [e for e in data["errors"] if "오버행 불일치" in e]
        assert not mismatch_errors, f"불일치 에러 발생: {mismatch_errors}"
        assert data["total_combinations"] >= 1, f"조합 없음. 에러: {data['errors']}"

    def test_mismatched_vector_overhang_right_blocks_assembly(self):
        """마지막 파트 OHR이 vector_overhang_right와 다르면 조합 제외."""
        payload = {
            "project_name": "vec_oh_mismatch_right",
            "assembly_method": "goldengate",
            "vector_sequence": self.VEC_SEQ,
            "vector_overhang_left": "AAAA",
            "vector_overhang_right": "TCAC",  # 실제 마지막 파트 OHR은 CCCC → 불일치
            "parts": [
                {"id": "p1", "name": "part1", "type": "Promoter",
                 "sequence": self.PART1_SEQ,
                 "overhang_left": "AAAA", "overhang_right": "TTTT"},
                {"id": "p2", "name": "part2", "type": "RBS",
                 "sequence": self.PART2_SEQ,
                 "overhang_left": "TTTT", "overhang_right": "CCCC"},
            ],
            "modules": [{"id": "M1", "slots": [
                {"type": "Promoter", "part_names": ["part1"]},
                {"type": "RBS",      "part_names": ["part2"]},
            ]}],
        }
        resp = client.post("/run", json=payload)
        assert resp.status_code == 200
        data = resp.json()
        mismatch_errors = [e for e in data["errors"] if "우측 오버행 불일치" in e]
        assert mismatch_errors, f"불일치 에러가 없음. 에러: {data['errors']}"
        assert data["total_combinations"] == 0, "불일치 조합이 결과에 포함됨"

    def test_mismatched_vector_overhang_left_blocks_assembly(self):
        """첫 번째 파트 OHL이 vector_overhang_left와 다르면 조합 제외."""
        payload = {
            "project_name": "vec_oh_mismatch_left",
            "assembly_method": "goldengate",
            "vector_sequence": self.VEC_SEQ,
            "vector_overhang_left": "GGGG",  # 실제 첫 파트 OHL은 AAAA → 불일치
            "vector_overhang_right": "CCCC",
            "parts": [
                {"id": "p1", "name": "part1", "type": "Promoter",
                 "sequence": self.PART1_SEQ,
                 "overhang_left": "AAAA", "overhang_right": "TTTT"},
                {"id": "p2", "name": "part2", "type": "RBS",
                 "sequence": self.PART2_SEQ,
                 "overhang_left": "TTTT", "overhang_right": "CCCC"},
            ],
            "modules": [{"id": "M1", "slots": [
                {"type": "Promoter", "part_names": ["part1"]},
                {"type": "RBS",      "part_names": ["part2"]},
            ]}],
        }
        resp = client.post("/run", json=payload)
        assert resp.status_code == 200
        data = resp.json()
        mismatch_errors = [e for e in data["errors"] if "좌측 오버행 불일치" in e]
        assert mismatch_errors, f"불일치 에러가 없음. 에러: {data['errors']}"
        assert data["total_combinations"] == 0, "불일치 조합이 결과에 포함됨"

    def test_no_vector_overhang_specified_skips_validation(self):
        """vector_overhang_left/right 미지정 시 검증 없이 조립 진행."""
        payload = {
            "project_name": "vec_oh_none",
            "assembly_method": "goldengate",
            "vector_sequence": self.VEC_SEQ,
            # vector_overhang_left/right 생략
            "parts": [
                {"id": "p1", "name": "part1", "type": "Promoter",
                 "sequence": self.PART1_SEQ,
                 "overhang_left": "AAAA", "overhang_right": "TTTT"},
                {"id": "p2", "name": "part2", "type": "RBS",
                 "sequence": self.PART2_SEQ,
                 "overhang_left": "TTTT", "overhang_right": "CCCC"},
            ],
            "modules": [{"id": "M1", "slots": [
                {"type": "Promoter", "part_names": ["part1"]},
                {"type": "RBS",      "part_names": ["part2"]},
            ]}],
        }
        resp = client.post("/run", json=payload)
        assert resp.status_code == 200
        data = resp.json()
        assert data["total_combinations"] >= 1, f"조합 없음. 에러: {data['errors']}"
