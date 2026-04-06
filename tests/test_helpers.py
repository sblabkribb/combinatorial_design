"""
Unit tests for combinatorial design helper functions.
These test pure-logic functions that don't require a running server.
"""

import os
import sys
import tempfile
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

from api import (
    get_positions_by_label,
    remove_duplicate_features,
    record_to_genbank_b64,
    extract_features,
    create_insert_genbank_files,
    filter_combinations_by_overhang,
    PartInput,
    ModuleInput,
    ModuleSlot,
)


# ── Fixtures ──────────────────────────────────────────────────────────────────

def make_record_with_labels(labels_and_positions: list[tuple]) -> SeqRecord:
    """Create a SeqRecord with misc_features at given (label, start, end) positions."""
    seq_len = max(end for _, _, end in labels_and_positions) + 10
    record = SeqRecord(Seq("A" * seq_len), id="test", name="test")
    record.annotations = {"molecule_type": "DNA"}
    for label, start, end in labels_and_positions:
        feature = SeqFeature(
            FeatureLocation(start, end, strand=1),
            type="misc_feature",
            qualifiers={"label": [label]},
        )
        record.features.append(feature)
    return record


# ── get_positions_by_label ────────────────────────────────────────────────────

class TestGetPositionsByLabel:
    def test_single_label_found(self):
        record = make_record_with_labels([("promoter", 10, 50)])
        result = get_positions_by_label(record, ["promoter"])
        assert result == {"start": 10, "end": 50}

    def test_multiple_labels_returns_bounding_box(self):
        record = make_record_with_labels([("O1", 5, 20), ("CDS", 20, 80), ("O2", 80, 95)])
        result = get_positions_by_label(record, ["O1", "CDS", "O2"])
        assert result == {"start": 5, "end": 95}

    def test_label_not_found_returns_none(self):
        record = make_record_with_labels([("CDS", 10, 50)])
        result = get_positions_by_label(record, ["promoter"])
        assert result == {"start": None, "end": None}

    def test_partial_match_returns_found_only(self):
        record = make_record_with_labels([("CDS", 10, 50)])
        result = get_positions_by_label(record, ["CDS", "BsaI"])
        # Only CDS found → start=10, end=50
        assert result == {"start": 10, "end": 50}

    def test_empty_record(self):
        record = SeqRecord(Seq("ATCG"), id="empty")
        result = get_positions_by_label(record, ["CDS"])
        assert result == {"start": None, "end": None}


# ── remove_duplicate_features ─────────────────────────────────────────────────

class TestRemoveDuplicateFeatures:
    def test_no_duplicates_unchanged(self):
        record = make_record_with_labels([("a", 0, 10), ("b", 10, 20)])
        result = remove_duplicate_features(record)
        assert len(result.features) == 2

    def test_duplicate_removed(self):
        record = make_record_with_labels([("a", 0, 10), ("a_dup", 0, 10)])
        result = remove_duplicate_features(record)
        assert len(result.features) == 1

    def test_first_occurrence_kept(self):
        record = make_record_with_labels([("original", 0, 10), ("copy", 0, 10)])
        result = remove_duplicate_features(record)
        label = result.features[0].qualifiers["label"][0]
        assert label == "original"


# ── record_to_genbank_b64 ─────────────────────────────────────────────────────

class TestRecordToGenbankB64:
    def test_returns_valid_base64(self):
        import base64
        record = SeqRecord(Seq("ATCGATCG"), id="test", name="test",
                           annotations={"molecule_type": "DNA"})
        b64 = record_to_genbank_b64(record)
        decoded = base64.b64decode(b64).decode()
        assert "LOCUS" in decoded
        assert "ATCGATCG" in decoded.upper() or "atcgatcg" in decoded.lower()

    def test_roundtrip_preserves_features(self):
        import base64
        from Bio import SeqIO
        import io
        record = make_record_with_labels([("test_feat", 0, 8)])
        record.seq = Seq("ATCGATCG")
        b64 = record_to_genbank_b64(record)
        decoded = base64.b64decode(b64).decode()
        loaded = SeqIO.read(io.StringIO(decoded), "genbank")
        assert len(loaded.features) == 1


# ── create_insert_genbank_files ───────────────────────────────────────────────

class TestCreateInsertGenbankFiles:
    def test_creates_file_per_part(self):
        parts = [
            PartInput(id="p1", name="T7prom", type="Promoter",
                      sequence="TAATACGACTCACTATA",
                      overhang_left="O1", overhang_right="O2"),
            PartInput(id="p2", name="RBS1", type="RBS",
                      sequence="AAAGAGGAGAAA",
                      overhang_left="O2", overhang_right="O3"),
        ]
        with tempfile.TemporaryDirectory() as tmpdir:
            create_insert_genbank_files(parts, tmpdir)
            files = os.listdir(tmpdir)
            assert "O1-T7prom-O2-insert.gb" in files
            assert "O2-RBS1-O3-insert.gb" in files

    def test_file_contains_correct_sequence(self):
        from Bio import SeqIO
        parts = [
            PartInput(id="p1", name="GFP", type="CDS",
                      sequence="ATGAGTAAAGGAGAAGAA",
                      overhang_left="O3", overhang_right="O4"),
        ]
        with tempfile.TemporaryDirectory() as tmpdir:
            create_insert_genbank_files(parts, tmpdir)
            record = SeqIO.read(
                os.path.join(tmpdir, "O3-GFP-O4-insert.gb"), "genbank"
            )
            assert str(record.seq).upper() == "ATGAGTAAAGGAGAAGAA"

    def test_reverse_direction_applies_revcomp(self):
        from Bio import SeqIO
        seq = "ATCGATCG"
        expected = str(Seq(seq).reverse_complement())
        parts = [
            PartInput(id="p1", name="Rev", type="CDS",
                      sequence=seq, direction="reverse",
                      overhang_left="O1", overhang_right="O2"),
        ]
        with tempfile.TemporaryDirectory() as tmpdir:
            create_insert_genbank_files(parts, tmpdir)
            record = SeqIO.read(
                os.path.join(tmpdir, "O1-Rev-O2-insert.gb"), "genbank"
            )
            assert str(record.seq).upper() == expected.upper()

    def test_feature_label_matches_part_name(self):
        from Bio import SeqIO
        parts = [
            PartInput(id="p1", name="BBa_J23100", type="Promoter",
                      sequence="TTGACAGCTAGCTCAGT",
                      overhang_left="O1", overhang_right="O2"),
        ]
        with tempfile.TemporaryDirectory() as tmpdir:
            create_insert_genbank_files(parts, tmpdir)
            record = SeqIO.read(
                os.path.join(tmpdir, "O1-BBa_J23100-O2-insert.gb"), "genbank"
            )
            labels = [f.qualifiers["label"][0] for f in record.features]
            assert "BBa_J23100" in labels


# ── filter_combinations_by_overhang ──────────────────────────────────────────

class TestFilterCombinationsByOverhang:
    def make_withvector_files(self, specs: list[tuple]) -> list[str]:
        """specs: [(overhang_left, name, overhang_right), ...]"""
        return [f"{ol}-{name}-{or_}-withvector.gb"
                for ol, name, or_ in specs]

    def test_compatible_combination_passes_filter(self):
        files = self.make_withvector_files([
            ("O1", "T7prom", "O2"),
            ("O2", "RBS1",   "O3"),
            ("O3", "GFP",    "O4"),
        ])
        modules = [ModuleInput(
            id="MM1",
            slots=[
                ModuleSlot(type="Promoter", part_names=["T7prom"]),
                ModuleSlot(type="RBS",      part_names=["RBS1"]),
                ModuleSlot(type="CDS",      part_names=["GFP"]),
            ],
        )]
        result = filter_combinations_by_overhang(modules, files)
        assert len(result) == 1

    def test_incompatible_overhang_filtered_out(self):
        # T7prom ends with O2, but RBS1 starts with O3 → mismatch
        files = self.make_withvector_files([
            ("O1", "T7prom", "O2"),
            ("O3", "RBS1",   "O4"),  # O3 != O2 → incompatible
        ])
        modules = [ModuleInput(
            id="MM1",
            slots=[
                ModuleSlot(type="Promoter", part_names=["T7prom"]),
                ModuleSlot(type="RBS",      part_names=["RBS1"]),
            ],
        )]
        result = filter_combinations_by_overhang(modules, files)
        assert len(result) == 0

    def test_combinatorial_expansion(self):
        # 2 promoters × 1 RBS → 2 combinations
        files = self.make_withvector_files([
            ("O1", "ProA", "O2"),
            ("O1", "ProB", "O2"),
            ("O2", "RBS1", "O3"),
        ])
        modules = [ModuleInput(
            id="MM1",
            slots=[
                ModuleSlot(type="Promoter", part_names=["ProA", "ProB"]),
                ModuleSlot(type="RBS",      part_names=["RBS1"]),
            ],
        )]
        result = filter_combinations_by_overhang(modules, files)
        assert len(result) == 2

    def test_three_way_expansion(self):
        # 2 promoters × 2 RBS × 1 CDS → 4 valid combinations
        files = self.make_withvector_files([
            ("O1", "ProA", "O2"), ("O1", "ProB", "O2"),
            ("O2", "RBS1", "O3"), ("O2", "RBS2", "O3"),
            ("O3", "GFP",  "O4"),
        ])
        modules = [ModuleInput(
            id="MM1",
            slots=[
                ModuleSlot(type="Promoter", part_names=["ProA", "ProB"]),
                ModuleSlot(type="RBS",      part_names=["RBS1", "RBS2"]),
                ModuleSlot(type="CDS",      part_names=["GFP"]),
            ],
        )]
        result = filter_combinations_by_overhang(modules, files)
        assert len(result) == 4
