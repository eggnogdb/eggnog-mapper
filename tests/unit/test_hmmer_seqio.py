"""Unit tests for eggnogmapper.search.hmmer.hmmer_seqio.iter_fasta_seqs.

Focus: comment lines before the first record. '#' lines were already
skipped; ';' lines (the original Pearson FASTA comment convention, also used
by FAIR-bioHeaders ';~' metadata headers) are skipped the same way. Any other
text before the first '>' must still be rejected.
"""

import bz2
import gzip

import pytest

from eggnogmapper.search.hmmer.hmmer_seqio import iter_fasta_seqs


RECORDS = (
    ">seq1 some description\n"
    "MKVLAAGIVG\n"
    "LLAWQ\n"
    ">seq2\n"
    "MSTNPKPQRK\n"
)
EXPECTED = [("seq1", "MKVLAAGIVGLLAWQ"), ("seq2", "MSTNPKPQRK")]


def _write(tmp_path, text, name="in.fa"):
    path = tmp_path / name
    path.write_text(text)
    return str(path)


def _parse(source, **kwargs):
    return list(iter_fasta_seqs(source, silent=True, **kwargs))


def test_plain_fasta(tmp_path):
    assert _parse(_write(tmp_path, RECORDS)) == EXPECTED


def test_hash_comment_before_first_record(tmp_path):
    assert _parse(_write(tmp_path, "# comment\n" + RECORDS)) == EXPECTED


def test_semicolon_comment_before_first_record(tmp_path):
    assert _parse(_write(tmp_path, ";comment\n" + RECORDS)) == EXPECTED


def test_fair_bioheaders_header_before_first_record(tmp_path):
    header = (
        ";~schema: https://w3id.org/fair-bioheaders/fhp/v0.4.0\n"
        ";~schemaVersion: '1.0'\n"
        ";~taxon:\n"
        ";~  name: Homo sapiens\n"
        ";~checksum: NJ3wlm9hfwvHy0PV2ytmd7Ey4El+lBCkQFnTEc3M2eE=\n"
    )
    assert _parse(_write(tmp_path, header + RECORDS)) == EXPECTED


def test_mixed_comments_and_blank_lines_before_first_record(tmp_path):
    text = ";one\n\n# two\n;\n" + RECORDS
    assert _parse(_write(tmp_path, text)) == EXPECTED


def test_semicolon_comment_crlf(tmp_path):
    path = tmp_path / "crlf.fa"
    path.write_bytes((";comment\n" + RECORDS).replace("\n", "\r\n").encode())
    assert _parse(str(path)) == EXPECTED


def test_semicolon_comment_gzip(tmp_path):
    path = tmp_path / "in.fa.gz"
    with gzip.open(path, "wt") as fh:
        fh.write(";~schema: x\n" + RECORDS)
    assert _parse(str(path)) == EXPECTED


def test_semicolon_comment_bz2(tmp_path):
    path = tmp_path / "in.fa.bz2"
    with bz2.open(path, "wt") as fh:
        fh.write(";~schema: x\n" + RECORDS)
    assert _parse(str(path)) == EXPECTED


def test_semicolon_comment_string_source():
    assert _parse(";comment\n" + RECORDS) == EXPECTED


def test_parsing_after_first_record_unchanged(tmp_path):
    # Only lines before the first record are skipped; later lines are
    # handled exactly as before this change.
    text = ";comment\n>a\nMKV\n;c\n>b\nMST\n"
    assert _parse(_write(tmp_path, text)) == [("a", "MKV;c"), ("b", "MST")]


def test_semicolon_comment_translate(tmp_path):
    text = ";comment\n>cds1\nATGAAAGTTTAA\n"
    assert _parse(_write(tmp_path, text), translate=True) == [("cds1", "MKV")]


@pytest.mark.parametrize("text", ["", ";comment\n# comment\n"])
def test_empty_or_comments_only_yields_nothing(tmp_path, text):
    assert _parse(_write(tmp_path, text)) == []


@pytest.mark.parametrize("prefix", ["", ";comment\n"])
@pytest.mark.parametrize("junk", ["MKVLAAGIVG", "random text", "@read1"])
def test_other_text_before_first_record_still_fails(tmp_path, prefix, junk):
    with pytest.raises(Exception, match="Wrong format"):
        _parse(_write(tmp_path, prefix + junk + "\n" + RECORDS))
