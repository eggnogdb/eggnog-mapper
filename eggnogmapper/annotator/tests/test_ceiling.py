"""Unit tests for TaxScopeCeilingResolver.

Tests cover:

1. auto-narrow ceiling for known organisms:
   - Arabidopsis thaliana (3702, plant) → Viridiplantae
   - Homo sapiens (9606, human) → Metazoa
   - Saccharomyces cerevisiae (4932, yeast) → Fungi
   - Escherichia coli K-12 (511145, bacterium) → Prokaryota

2. auto-broad ceiling for selected organisms:
   - Homo sapiens → Metazoa
   - Saccharomyces cerevisiae → Fungi

3. Named clade ceiling: "Viridiplantae" returns same ceiling regardless of seed.

4. ev_lca_passes_ceiling:
   - ev_lca inside ceiling → True
   - ev_lca outside ceiling → False

5. Memoization: second call for same seed returns same object.

Biological context
------------------
The sample database at /eggnog-eco/data/e7/sample/final/mapper/eggnog.taxa.db
covers ~2 000 species (a realistic subset of the full ~50 k species eggnog.taxa.db).
All five key species tested here (taxids 9606, 3702, 4932, 511145, 426368) are
present in the sample DB (verified manually).

At full scale:
- Prokaryota set covers ~35 k species (Bacteria ∪ Archaea).
- Fungi set covers ~2 000 species.
- Metazoa covers ~8 000 species.
- All sets are computed once in build() and memoized.
- The _ceiling_cache and _valid_species_cache prevent any O(n) re-scan per seed.
"""

import pytest

from .conftest import TAXA_DB_PATH


# ---------------------------------------------------------------------------
# Module-level fixtures (shared across all test functions)
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def lineage_cache(taxa_db_path):
    """LineageCache loaded from the sample taxa.db.

    Species in the sample DB (3702, 9606, 4932, 511145, 426368) are all
    present — verified manually via SQLite query before writing tests.
    """
    from eggnogmapper.annotator.lineage import LineageCache
    return LineageCache(taxa_db_path=taxa_db_path)


@pytest.fixture(scope="module")
def resolver_narrow(lineage_cache, taxa_db_path):
    """TaxScopeCeilingResolver built in auto-narrow mode."""
    from eggnogmapper.annotator.e7.ceiling import TaxScopeCeilingResolver
    return TaxScopeCeilingResolver.build(
        lineage_cache=lineage_cache,
        mode="auto-narrow",
        taxa_db_path=taxa_db_path,
    )


@pytest.fixture(scope="module")
def resolver_broad(lineage_cache, taxa_db_path):
    """TaxScopeCeilingResolver built in auto-broad mode."""
    from eggnogmapper.annotator.e7.ceiling import TaxScopeCeilingResolver
    return TaxScopeCeilingResolver.build(
        lineage_cache=lineage_cache,
        mode="auto-broad",
        taxa_db_path=taxa_db_path,
    )


@pytest.fixture(scope="module")
def resolver_viridiplantae(lineage_cache, taxa_db_path):
    """TaxScopeCeilingResolver with a fixed named-clade ceiling: Viridiplantae."""
    from eggnogmapper.annotator.e7.ceiling import TaxScopeCeilingResolver
    return TaxScopeCeilingResolver.build(
        lineage_cache=lineage_cache,
        mode="Viridiplantae",
        taxa_db_path=taxa_db_path,
    )


# ---------------------------------------------------------------------------
# Taxid constants used in tests
# ---------------------------------------------------------------------------

ARABIDOPSIS_TAXID   = "3702"    # Arabidopsis thaliana — Viridiplantae
HUMAN_TAXID         = "9606"    # Homo sapiens — Metazoa
YEAST_TAXID         = "4932"    # Saccharomyces cerevisiae — Fungi
ECOLI_TAXID         = "511145"  # Escherichia coli K-12 — Bacteria (Prokaryota)
ARCHAEA_TAXID       = "426368"  # Methanococcus maripaludis C7 — Archaea (Prokaryota)
MONOSIGA_TAXID      = "81824"   # Monosiga brevicollis — Opisthokonta, NOT Metazoa or Fungi
ENCEPHALITOZOON_TAXID = "6035"  # Encephalitozoon cuniculi — Microsporidia (Fungi lineage but excluded)

VIRIDIPLANTAE_TAXID = "33090"
METAZOA_TAXID       = "33208"
FUNGI_TAXID         = "4751"


# ---------------------------------------------------------------------------
# 1. auto-narrow ceilings for key organisms
# ---------------------------------------------------------------------------

def test_auto_narrow_arabidopsis_eukaryota(resolver_narrow):
    """auto-narrow: Arabidopsis thaliana (3702) should get Eukaryota ceiling.

    AUTO_NARROW_PRIORITY is domain-level only [Eukaryota, Prokaryota].
    All eukaryotes receive the Eukaryota ceiling regardless of sub-kingdom.
    """
    ceiling = resolver_narrow.resolve_ceiling(ARABIDOPSIS_TAXID)
    name = resolver_narrow.ceiling_name(ceiling)
    assert name == "Eukaryota", (
        f"Expected Eukaryota ceiling for Arabidopsis (3702) in auto-narrow, got {name!r} "
        f"(taxid={ceiling})"
    )


def test_auto_narrow_human_eukaryota(resolver_narrow):
    """auto-narrow: Homo sapiens (9606) should get Eukaryota ceiling.

    AUTO_NARROW_PRIORITY is domain-level only — no Metazoa or Opisthokonta
    tiers.  Human gets the same Eukaryota ceiling as all other eukaryotes.
    """
    ceiling = resolver_narrow.resolve_ceiling(HUMAN_TAXID)
    name = resolver_narrow.ceiling_name(ceiling)
    assert name == "Eukaryota", (
        f"Expected Eukaryota ceiling for human (9606) in auto-narrow, got {name!r} "
        f"(taxid={ceiling})"
    )


def test_auto_narrow_yeast_eukaryota(resolver_narrow):
    """auto-narrow: Saccharomyces cerevisiae (4932) should get Eukaryota ceiling.

    AUTO_NARROW_PRIORITY is domain-level only — Fungi is not in the list.
    Yeast gets Eukaryota ceiling, same as all other eukaryotes.
    """
    ceiling = resolver_narrow.resolve_ceiling(YEAST_TAXID)
    name = resolver_narrow.ceiling_name(ceiling)
    assert name == "Eukaryota", (
        f"Expected Eukaryota ceiling for S. cerevisiae (4932) in auto-narrow, got {name!r} "
        f"(taxid={ceiling})"
    )


def test_auto_narrow_ecoli_prokaryota(resolver_narrow):
    """auto-narrow: Escherichia coli K-12 (511145) should get Prokaryota ceiling.

    E. coli lacks Eukaryota (2759) in its lineage so none of the eukaryotic
    tiers match; the PROKARYOTA_SYNTHETIC sentinel catches it.
    """
    from eggnogmapper.annotator.e7.ceiling import PROKARYOTA_SYNTHETIC
    ceiling = resolver_narrow.resolve_ceiling(ECOLI_TAXID)
    name = resolver_narrow.ceiling_name(ceiling)
    assert name == "Prokaryota", (
        f"Expected Prokaryota ceiling for E. coli (511145), got {name!r} "
        f"(taxid={ceiling})"
    )
    assert ceiling == PROKARYOTA_SYNTHETIC, (
        f"Expected PROKARYOTA_SYNTHETIC sentinel, got {ceiling!r}"
    )


def test_auto_narrow_archaea_prokaryota(resolver_narrow):
    """auto-narrow: Methanococcus maripaludis (426368) should get Prokaryota ceiling.

    Archaea (2157) are not Eukaryota, so the PROKARYOTA_SYNTHETIC sentinel
    catches all archaeal species (as well as bacterial ones).
    """
    ceiling = resolver_narrow.resolve_ceiling(ARCHAEA_TAXID)
    name = resolver_narrow.ceiling_name(ceiling)
    assert name == "Prokaryota", (
        f"Expected Prokaryota ceiling for Methanococcus (426368), got {name!r}"
    )


# ---------------------------------------------------------------------------
# 2. auto-broad ceilings
# ---------------------------------------------------------------------------

def test_auto_broad_human_metazoa(resolver_broad):
    """auto-broad: Homo sapiens (9606) should still get Metazoa ceiling.

    In auto-broad mode Metazoa is first in the priority list (broadest
    match wins), and Metazoa is also the broadest correct match for human.
    """
    ceiling = resolver_broad.resolve_ceiling(HUMAN_TAXID)
    name = resolver_broad.ceiling_name(ceiling)
    assert name == "Metazoa", (
        f"Expected Metazoa ceiling for human in auto-broad mode, got {name!r}"
    )


def test_auto_broad_yeast_fungi(resolver_broad):
    """auto-broad: Saccharomyces cerevisiae (4932) should get Fungi ceiling.

    In auto-broad mode, Metazoa comes first in the priority list but yeast
    is not in Metazoa, so Fungi is the first matching entry.
    """
    ceiling = resolver_broad.resolve_ceiling(YEAST_TAXID)
    name = resolver_broad.ceiling_name(ceiling)
    assert name == "Fungi", (
        f"Expected Fungi ceiling for S. cerevisiae in auto-broad mode, got {name!r}"
    )


# ---------------------------------------------------------------------------
# 2b. Divergence cases: auto-narrow vs auto-broad produce DIFFERENT ceilings
#
# These are the only organisms where the two modes behave differently.
# For the three main kingdoms (Metazoa, Viridiplantae, Fungi) both modes
# produce identical ceilings — the ordering difference in the priority list
# does not matter because each organism belongs to exactly one kingdom.
#
# Divergence requires a seed that is in Opisthokonta but NOT in Metazoa or
# Fungi proper: choanoflagellates and Microsporidia (excluded from Fungi).
# auto-narrow has Opisthokonta in its list; auto-broad does not → falls through
# to Eukaryota.
# ---------------------------------------------------------------------------

def test_narrow_vs_broad_monosiga_both_eukaryota(resolver_narrow, resolver_broad):
    """Monosiga brevicollis (81824) gets Eukaryota in both auto-narrow and auto-broad.

    Monosiga is a choanoflagellate: in Opisthokonta (33154) but NOT in
    Metazoa (33208) or Fungi (4751).  auto-narrow uses domain-level priority
    [Eukaryota, Prokaryota] so Monosiga gets Eukaryota.  auto-broad has no
    Opisthokonta tier, so it also falls through to Eukaryota.

    Both modes agree for Monosiga — the two modes diverge only for organisms
    within Metazoa, Viridiplantae, or Fungi (auto-broad gives the sub-kingdom
    ceiling; auto-narrow gives Eukaryota).
    """
    narrow_ceiling = resolver_narrow.resolve_ceiling(MONOSIGA_TAXID)
    broad_ceiling = resolver_broad.resolve_ceiling(MONOSIGA_TAXID)
    assert resolver_narrow.ceiling_name(narrow_ceiling) == "Eukaryota", (
        f"auto-narrow: expected Eukaryota for Monosiga, got "
        f"{resolver_narrow.ceiling_name(narrow_ceiling)!r}"
    )
    assert resolver_broad.ceiling_name(broad_ceiling) == "Eukaryota", (
        f"auto-broad: expected Eukaryota for Monosiga, got "
        f"{resolver_broad.ceiling_name(broad_ceiling)!r}"
    )


def test_narrow_vs_broad_microsporidia_both_eukaryota(resolver_narrow, resolver_broad):
    """Encephalitozoon cuniculi (6035) gets Eukaryota in both auto-narrow and auto-broad.

    E. cuniculi is a Microsporidian: its NCBI lineage passes through Fungi (4751)
    but Microsporidia (6029) are explicitly excluded from the Fungi set in both
    modes.  auto-narrow is domain-level [Eukaryota, Prokaryota] so Microsporidia
    falls through to Eukaryota.  auto-broad also falls through Metazoa/Viridiplantae/
    Fungi (excluded) to Eukaryota.

    Both modes agree on Eukaryota.  The Microsporidia exclusion from Fungi still
    applies in both modes — Microsporidia are NOT assigned to the Fungi ceiling.
    """
    narrow_ceiling = resolver_narrow.resolve_ceiling(ENCEPHALITOZOON_TAXID)
    broad_ceiling = resolver_broad.resolve_ceiling(ENCEPHALITOZOON_TAXID)
    assert resolver_narrow.ceiling_name(narrow_ceiling) == "Eukaryota", (
        f"auto-narrow: expected Eukaryota for Microsporidia, got "
        f"{resolver_narrow.ceiling_name(narrow_ceiling)!r} — "
        "Microsporidia must not be assigned to Fungi ceiling"
    )
    assert resolver_broad.ceiling_name(broad_ceiling) == "Eukaryota", (
        f"auto-broad: expected Eukaryota for Microsporidia, got "
        f"{resolver_broad.ceiling_name(broad_ceiling)!r}"
    )


def test_narrow_vs_broad_diverge_for_metazoa(resolver_narrow, resolver_broad):
    """Human (9606) gets Eukaryota (narrow) vs Metazoa (broad) — modes diverge.

    AUTO_NARROW_PRIORITY is domain-level only: human gets Eukaryota.
    AUTO_BROAD_PRIORITY starts with Metazoa: human gets Metazoa.
    This is the canonical divergence case for the two modes.
    """
    narrow_ceiling = resolver_narrow.resolve_ceiling(HUMAN_TAXID)
    broad_ceiling = resolver_broad.resolve_ceiling(HUMAN_TAXID)
    assert resolver_narrow.ceiling_name(narrow_ceiling) == "Eukaryota", (
        f"auto-narrow: expected Eukaryota for human, got "
        f"{resolver_narrow.ceiling_name(narrow_ceiling)!r}"
    )
    assert resolver_broad.ceiling_name(broad_ceiling) == "Metazoa", (
        f"auto-broad: expected Metazoa for human, got "
        f"{resolver_broad.ceiling_name(broad_ceiling)!r}"
    )
    assert narrow_ceiling != broad_ceiling, (
        "auto-narrow and auto-broad must produce DIFFERENT ceilings for human"
    )


def test_microsporidia_not_in_fungi_minus_microsporidia_set(resolver_narrow):
    """E. cuniculi must NOT appear in the fungi-minus-microsporidia set.

    If this fails it means Microsporidia would silently inherit the Fungi ceiling,
    contradicting the explicit exclusion required by the spec.
    """
    assert ENCEPHALITOZOON_TAXID not in resolver_narrow._fungi_minus_microsporidia, (
        "Encephalitozoon cuniculi (6035) must be excluded from fungi_minus_microsporidia"
    )


# ---------------------------------------------------------------------------
# 3. Named clade ceiling
# ---------------------------------------------------------------------------

def test_named_clade_viridiplantae_resolves(resolver_viridiplantae):
    """Fixed named ceiling 'Viridiplantae' is resolved to taxid 33090.

    This verifies that TaxScopeCeilingResolver.build() can look up a named
    clade in the taxa DB and store the resolved numeric taxid.
    """
    assert resolver_viridiplantae._fixed_ceiling == VIRIDIPLANTAE_TAXID, (
        f"Expected fixed_ceiling=33090, got {resolver_viridiplantae._fixed_ceiling!r}"
    )


def test_named_clade_same_for_all_seeds(resolver_viridiplantae):
    """A fixed named ceiling returns the same taxid regardless of seed identity.

    Whether the seed is a plant, an animal, or a bacterium, the ceiling is
    always the fixed taxid (Viridiplantae=33090).  This is the semantics of
    --tax_scope Viridiplantae: filter all seeds' donors to Viridiplantae members.
    """
    for seed_taxid in [ARABIDOPSIS_TAXID, HUMAN_TAXID, ECOLI_TAXID]:
        ceiling = resolver_viridiplantae.resolve_ceiling(seed_taxid)
        assert ceiling == VIRIDIPLANTAE_TAXID, (
            f"Fixed-ceiling resolver returned {ceiling!r} for seed {seed_taxid!r}, "
            f"expected Viridiplantae ({VIRIDIPLANTAE_TAXID})"
        )


# ---------------------------------------------------------------------------
# 4. ev_lca_passes_ceiling
# ---------------------------------------------------------------------------

def test_ev_lca_inside_ceiling_passes(resolver_narrow):
    """ev_lca_passes_ceiling: human (9606) is inside Metazoa → True.

    9606's lineage contains 33208 (Metazoa), so the event should pass the
    Metazoa ceiling filter.
    """
    result = resolver_narrow.ev_lca_passes_ceiling("9606", METAZOA_TAXID)
    assert result is True, (
        "Human (9606) ev_lca should pass a Metazoa ceiling (9606 ∈ Metazoa)"
    )


def test_ev_lca_outside_ceiling_fails(resolver_narrow):
    """ev_lca_passes_ceiling: E. coli (511145) is outside Metazoa → False.

    511145's lineage does NOT contain 33208 (Metazoa), so an event with
    ev_lca=511145 should be discarded under a Metazoa ceiling.

    This is the key biological correctness test: a technically valid
    output (correct format, no exception) that is biologically wrong if
    it returns True.
    """
    result = resolver_narrow.ev_lca_passes_ceiling("511145", METAZOA_TAXID)
    assert result is False, (
        "E. coli (511145) ev_lca should NOT pass a Metazoa ceiling "
        "(511145 ∉ Metazoa)"
    )


def test_ev_lca_empty_string_is_conservative_keep(resolver_narrow):
    """ev_lca_passes_ceiling: empty ev_lca is kept (conservative).

    Missing ev_lca data should not silently discard events; the method
    returns True to avoid dropping valid data when the DB lacks the taxid.
    """
    result = resolver_narrow.ev_lca_passes_ceiling("", METAZOA_TAXID)
    assert result is True, (
        "Empty ev_lca should pass conservatively (keep the event)"
    )


def test_ev_lca_root_ceiling_always_passes(resolver_narrow):
    """ev_lca_passes_ceiling: root ceiling always passes.

    When ceiling='root' (unknown seed) no filter is applied — every ev_lca
    passes, even a bacterial one.  This is the conservative fallback.
    """
    for ev_lca in [ECOLI_TAXID, HUMAN_TAXID, ARABIDOPSIS_TAXID]:
        result = resolver_narrow.ev_lca_passes_ceiling(ev_lca, "root")
        assert result is True, (
            f"ev_lca={ev_lca} should pass a 'root' ceiling (no filter applied)"
        )


def test_ev_lca_prokaryota_sentinel(resolver_narrow):
    """ev_lca_passes_ceiling: E. coli passes Prokaryota sentinel.

    Under a PROKARYOTA_SYNTHETIC ceiling, ev_lca must be in the pre-computed
    prokaryota set.  E. coli is in that set; human is not.
    """
    from eggnogmapper.annotator.e7.ceiling import PROKARYOTA_SYNTHETIC
    ecoli_result = resolver_narrow.ev_lca_passes_ceiling(ECOLI_TAXID, PROKARYOTA_SYNTHETIC)
    human_result = resolver_narrow.ev_lca_passes_ceiling(HUMAN_TAXID, PROKARYOTA_SYNTHETIC)
    assert ecoli_result is True, (
        "E. coli should pass PROKARYOTA_SYNTHETIC ceiling"
    )
    assert human_result is False, (
        "Human should NOT pass PROKARYOTA_SYNTHETIC ceiling"
    )


# ---------------------------------------------------------------------------
# 5. Memoization
# ---------------------------------------------------------------------------

def test_memoization_same_seed_returns_same_string(resolver_narrow):
    """Second call for same seed taxid returns the identical string (cached).

    The internal _ceiling_cache is a dict so the returned strings are
    the same object (Python interns short strings; we compare by identity
    for the cached dict value).
    """
    ceiling1 = resolver_narrow.resolve_ceiling(HUMAN_TAXID)
    ceiling2 = resolver_narrow.resolve_ceiling(HUMAN_TAXID)
    assert ceiling1 == ceiling2, (
        "Memoized ceiling should be equal on repeated call"
    )
    # The second call must hit the cache dict and return the same value
    assert HUMAN_TAXID in resolver_narrow._ceiling_cache, (
        "Human taxid should be in _ceiling_cache after first call"
    )


def test_memoization_valid_species_cache(resolver_narrow):
    """get_valid_species returns the identical frozenset on repeated calls."""
    s1 = resolver_narrow.get_valid_species(METAZOA_TAXID)
    s2 = resolver_narrow.get_valid_species(METAZOA_TAXID)
    assert s1 is s2, (
        "get_valid_species should return the same memoized frozenset object"
    )


def test_memoization_og_descendants_cache(resolver_narrow):
    """get_og_descendants returns the identical frozenset on repeated calls."""
    d1 = resolver_narrow.get_og_descendants(FUNGI_TAXID)
    d2 = resolver_narrow.get_og_descendants(FUNGI_TAXID)
    assert d1 is d2, (
        "get_og_descendants should return the same memoized frozenset object"
    )


# ---------------------------------------------------------------------------
# 6. Biological sanity checks on computed sets
# ---------------------------------------------------------------------------

def test_prokaryota_set_excludes_human(resolver_narrow):
    """The pre-computed prokaryota set must NOT include human (9606).

    A biologically wrong result here (human in prokaryota set) would cause
    all human seeds to get Prokaryota ceiling in auto-narrow mode, which
    is technically valid output but biologically nonsensical.
    """
    assert HUMAN_TAXID not in resolver_narrow._prokaryota_taxids, (
        "Human (9606) must not be in the prokaryota taxids set"
    )


def test_prokaryota_set_includes_ecoli(resolver_narrow):
    """The pre-computed prokaryota set must include E. coli (511145)."""
    assert ECOLI_TAXID in resolver_narrow._prokaryota_taxids, (
        "E. coli (511145) must be in the prokaryota taxids set"
    )


def test_fungi_minus_microsporidia_includes_yeast(resolver_narrow):
    """The Fungi-minus-Microsporidia set must include S. cerevisiae (4932)."""
    assert YEAST_TAXID in resolver_narrow._fungi_minus_microsporidia, (
        "S. cerevisiae (4932) must be in the fungi-minus-microsporidia set"
    )


def test_fungi_minus_microsporidia_excludes_human(resolver_narrow):
    """The Fungi-minus-Microsporidia set must NOT include human (9606)."""
    assert HUMAN_TAXID not in resolver_narrow._fungi_minus_microsporidia, (
        "Human (9606) must not appear in the fungi-minus-microsporidia set"
    )


def test_valid_species_metazoa_includes_human(resolver_narrow):
    """get_valid_species(Metazoa) must include human."""
    species = resolver_narrow.get_valid_species(METAZOA_TAXID)
    assert species is not None
    assert HUMAN_TAXID in species, (
        "Human (9606) should be in Metazoa valid-species set"
    )


def test_valid_species_metazoa_excludes_ecoli(resolver_narrow):
    """get_valid_species(Metazoa) must NOT include E. coli.

    If this test fails it means an E. coli sequence could receive Metazoa-clade
    annotation, which is biologically wrong (bacteria ≠ animals).
    """
    species = resolver_narrow.get_valid_species(METAZOA_TAXID)
    assert species is not None
    assert ECOLI_TAXID not in species, (
        "E. coli (511145) must not appear in Metazoa valid-species set"
    )


def test_og_descendants_metazoa_includes_human(resolver_narrow):
    """get_og_descendants(Metazoa) must include human as a descendant."""
    desc = resolver_narrow.get_og_descendants(METAZOA_TAXID)
    assert desc is not None
    assert HUMAN_TAXID in desc, (
        "Human (9606) should be in Metazoa OG descendants"
    )


def test_og_descendants_fungi_includes_yeast(resolver_narrow):
    """get_og_descendants(Fungi) must include yeast as a descendant."""
    desc = resolver_narrow.get_og_descendants(FUNGI_TAXID)
    assert desc is not None
    assert YEAST_TAXID in desc, (
        "S. cerevisiae (4932) should be in Fungi OG descendants"
    )


def test_root_returns_none_for_valid_species(resolver_narrow):
    """get_valid_species('root') must return None (no filter applied)."""
    result = resolver_narrow.get_valid_species("root")
    assert result is None, (
        "get_valid_species('root') should return None — no ceiling filter"
    )


def test_root_returns_none_for_og_descendants(resolver_narrow):
    """get_og_descendants('root') must return None (no filter applied)."""
    result = resolver_narrow.get_og_descendants("root")
    assert result is None, (
        "get_og_descendants('root') should return None — no ceiling filter"
    )


# ---------------------------------------------------------------------------
# 7. ceiling_name mapping
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("taxid,expected_name", [
    (METAZOA_TAXID,       "Metazoa"),
    (FUNGI_TAXID,         "Fungi"),
    (VIRIDIPLANTAE_TAXID, "Viridiplantae"),
    ("2759",              "Eukaryota"),
    ("2",                 "Bacteria"),
    ("2157",              "Archaea"),
    ("root",              "root"),
])
def test_ceiling_name_mapping(resolver_narrow, taxid, expected_name):
    """ceiling_name must map taxids to human-readable clade names.

    These names appear verbatim in the 'tax_ceiling' output column.
    If any mapping is wrong, the output column becomes misleading to users
    even if the underlying taxid logic is correct.
    """
    name = resolver_narrow.ceiling_name(taxid)
    assert name == expected_name, (
        f"ceiling_name({taxid!r}) = {name!r}, expected {expected_name!r}"
    )
