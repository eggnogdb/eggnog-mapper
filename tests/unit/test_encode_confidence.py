"""Unit tests for the compact annotation_confidence encoder and legend.

Covers the robustness invariants:
  RI-1  Field order = ANNOTATIONS_HEADER (single source of truth).
  RI-2  Tier chars derived via tier_name[0].lower() from TIER_CONFIDENCE.
  RI-4  None/empty conf_dict → fixed-width all-hyphen output.
  RI-5  Unknown conf_dict key → WARNING, never exception.
"""

import io
import logging

import pytest

from eggnogmapper.annotation.output import (
    ANNOTATIONS_HEADER,
    _TIER_LEGEND_ENTRIES,
    _legend_lines,
    encode_confidence,
    output_annotations_header,
)


@pytest.fixture
def full_conf_dict():
    """conf_dict covering every ANNOTATIONS_HEADER field with mixed tiers."""
    return {
        "Preferred_name": "high",
        "GOs": "low",
        "EC": "medium",
        "KEGG_ko": "high",
        "KEGG_Pathway": "low",
        "KEGG_Module": "medium",
        "KEGG_Reaction": "high",
        "KEGG_rclass": "low",
        "BRITE": "high",
        "KEGG_TC": "medium",
        "CAZy": "low",
        "BiGG_Reaction": "high",
        "PFAMs": "high",
    }


# ---------------------------------------------------------------------------
# T-1  shape + positions
# ---------------------------------------------------------------------------
def test_encode_confidence_shape_and_positions(full_conf_dict):
    result = encode_confidence(full_conf_dict, ANNOTATIONS_HEADER)
    assert len(result) == len(ANNOTATIONS_HEADER) == 13
    assert set(result) <= set("hml-")
    for field, expected_char in [
        ("Preferred_name", "h"),
        ("GOs", "l"),
        ("EC", "m"),
        ("PFAMs", "h"),
    ]:
        pos = ANNOTATIONS_HEADER.index(field)
        assert result[pos] == expected_char, (
            f"{field}@{pos}: got {result[pos]!r}, expected {expected_char!r}"
        )


# ---------------------------------------------------------------------------
# T-2  REORDER lockstep — the robustness ask
# ---------------------------------------------------------------------------
def test_encode_confidence_reorder_invariant(full_conf_dict):
    canonical = encode_confidence(full_conf_dict, ANNOTATIONS_HEADER)
    reordered_fields = list(reversed(ANNOTATIONS_HEADER))
    reordered = encode_confidence(full_conf_dict, reordered_fields)
    assert len(reordered) == len(reordered_fields)
    for i, field in enumerate(reordered_fields):
        canonical_pos = ANNOTATIONS_HEADER.index(field)
        assert reordered[i] == canonical[canonical_pos], (
            f"field {field!r}: reordered[{i}]={reordered[i]!r} "
            f"!= canonical[{canonical_pos}]={canonical[canonical_pos]!r}"
        )


# ---------------------------------------------------------------------------
# T-3  unknown-key WARNING via caplog
# ---------------------------------------------------------------------------
def test_encode_confidence_unknown_key_warning(caplog):
    bad_dict = {"Preferred_name": "high", "UNKNOWN_FIELD_XYZ": "medium"}
    with caplog.at_level(logging.WARNING, logger="eggnogmapper.annotation.output"):
        result = encode_confidence(bad_dict, ANNOTATIONS_HEADER)
    assert any("UNKNOWN_FIELD_XYZ" in rec.message for rec in caplog.records), (
        "expected a WARNING mentioning the unknown key"
    )
    assert len(result) == 13
    assert result[ANNOTATIONS_HEADER.index("Preferred_name")] == "h"


# ---------------------------------------------------------------------------
# T-4  None / empty → 13 hyphens
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("conf_input", [None, {}])
def test_encode_confidence_none_or_empty_is_all_hyphens(conf_input):
    result = encode_confidence(conf_input, ANNOTATIONS_HEADER)
    assert result == "-" * len(ANNOTATIONS_HEADER) == "-" * 13


# ---------------------------------------------------------------------------
# T-5  legend present/absent under --no_file_comments
# ---------------------------------------------------------------------------
def test_legend_present_when_comments_enabled():
    buf = io.StringIO()
    output_annotations_header(
        buf, no_file_comments=False, md5_field=False, print_header=True
    )
    content = buf.getvalue()
    assert "## annotation_confidence:" in content
    assert "## confidence codes:" in content
    assert "## confidence field order:" in content


def test_legend_absent_when_no_file_comments():
    buf = io.StringIO()
    output_annotations_header(
        buf, no_file_comments=True, md5_field=False, print_header=True
    )
    content = buf.getvalue()
    assert "## annotation_confidence:" not in content
    assert "## confidence codes:" not in content
    assert "## confidence field order:" not in content


# ---------------------------------------------------------------------------
# T-6  legend field-list == " ".join(ANNOTATIONS_HEADER)
# ---------------------------------------------------------------------------
def test_legend_field_order_matches_annotations_header():
    legend_entries = _TIER_LEGEND_ENTRIES + (("-", "not annotated"),)
    lines = _legend_lines(ANNOTATIONS_HEADER, legend_entries)
    assert len(lines) == 3
    field_line = lines[2]
    assert field_line == (
        "## confidence field order: " + " ".join(ANNOTATIONS_HEADER)
    )


def test_legend_codes_derived_from_tier_confidence():
    """Codes line must be derived from TIER_CONFIDENCE, not hardcoded."""
    from eggnogmapper.annotator.e7.constants import TIER_CONFIDENCE

    legend_entries = _TIER_LEGEND_ENTRIES + (("-", "not annotated"),)
    lines = _legend_lines(ANNOTATIONS_HEADER, legend_entries)
    codes_line = lines[1]
    for tier_name in TIER_CONFIDENCE.values():
        code = tier_name[0].lower()
        assert f"{code}={tier_name}" in codes_line, (
            f"legend missing {code}={tier_name}"
        )
    assert "-=not annotated" in codes_line
