"""Tests for chem-ops. Reference values are literature or RDKit-stable values for known drugs.

Run:  python -m pytest -q tests
"""

import sys
from pathlib import Path

import pytest
from rdkit import Chem

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "server"))

from chemops import alerts, draw, library, mol, props, similarity

EX = ROOT / "examples"
LIB = str(EX / "screening-hits.csv")

ASPIRIN = "CC(=O)Oc1ccccc1C(=O)O"
CAFFEINE = "Cn1c(=O)c2c(ncn2C)n(C)c1=O"
IMATINIB = "Cc1ccc(cc1Nc1nccc(n1)c1cccnc1)NC(=O)c1ccc(cc1)CN1CCN(C)CC1"
GEFITINIB = "COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1"
RHODANINE = "O=C1NC(=S)SC1=Cc1ccccc1"
TAXOL = ("CC1=C2C(C(=O)C3(C(CC4C(C3C(C(C2(C)C)(CC1OC(=O)C(C(c1ccccc1)NC(=O)c1ccccc1)O)O)"
         "OC(=O)c1ccccc1)(CO4)OC(=O)C)O)C)OC(=O)C")


# ------------------------------------------------------------------ parsing and standardization
def test_parse_and_canonicalize():
    assert mol.canonical(mol.parse(ASPIRIN)) == mol.canonical(mol.parse("OC(=O)c1ccccc1OC(C)=O"))
    assert mol.parse("this_is_not_a_smiles") is None
    assert mol.parse("") is None


def test_salt_is_stripped_to_parent():
    salt = mol.parse("CN(C)CCOC(c1ccccc1)c1ccccc1.Cl")
    parent, changes = mol.standardize(salt)
    assert "." not in mol.canonical(parent)
    assert any("largest fragment" in c for c in changes)


def test_standardized_mol_is_usable_for_fingerprints():
    # Regression: the standardizer returns molecules without ring info initialized.
    parent, _ = mol.standardize(mol.parse(GEFITINIB))
    assert similarity.fingerprints([parent])[0].GetNumOnBits() > 0
    assert props.descriptors(parent)["aromatic_rings"] == 3


def test_inchikey_matches_for_same_molecule_written_differently():
    a, _ = mol.standardize(mol.parse(ASPIRIN))
    b, _ = mol.standardize(mol.parse("OC(=O)c1ccccc1OC(C)=O"))
    assert mol.inchikey(a) == mol.inchikey(b) == "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"


# ------------------------------------------------------------------ descriptors and rules
@pytest.mark.parametrize("smiles,expected_mw,expected_formula", [
    (ASPIRIN, 180.16, "C9H8O4"),
    (CAFFEINE, 194.19, "C8H10N4O2"),
    (IMATINIB, 493.62, "C29H31N7O"),
])
def test_descriptors_match_reference_values(smiles, expected_mw, expected_formula):
    d = props.descriptors(mol.parse(smiles))
    assert d["mw"] == pytest.approx(expected_mw, abs=0.05)
    assert d["formula"] == expected_formula
    assert 0.0 <= d["qed"] <= 1.0


def test_aspirin_descriptor_detail():
    d = props.descriptors(mol.parse(ASPIRIN))
    assert d["hbd"] == 1
    assert d["tpsa"] == pytest.approx(63.6, abs=0.1)
    assert d["aromatic_rings"] == 1
    assert d["formal_charge"] == 0


def test_taxol_breaks_lipinski_on_mw_and_hba():
    m = mol.parse(TAXOL)
    d = props.descriptors(m)
    assert d["mw"] == pytest.approx(853.9, abs=0.2)
    ro5 = next(r for r in props.rule_sets(m, d) if r["rule"].startswith("Lipinski"))
    assert not ro5["pass"] and ro5["violations"] == 2
    assert any("MW" in c for c in ro5["failed_criteria"])


def test_lipinski_tolerates_one_violation():
    m = mol.parse(IMATINIB)
    ro5 = next(r for r in props.rule_sets(m) if r["rule"].startswith("Lipinski"))
    assert ro5["pass"]  # imatinib sits just inside, with cLogP near the limit


def test_fragment_passes_rule_of_three_and_fails_ghose():
    m = mol.parse("c1ccc2[nH]ccc2c1")  # indole
    rules = {r["rule"]: r for r in props.rule_sets(m)}
    assert rules["Rule of 3 (fragments)"]["pass"]
    assert not rules["Ghose"]["pass"]  # Ghose excludes small fragments by construction


@pytest.mark.skipif(not props.sa_available(), reason="RDKit Contrib SA_Score not available")
def test_synthetic_accessibility_ranks_taxol_above_aspirin():
    assert (props.descriptors(mol.parse(TAXOL))["sa_score"]
            > props.descriptors(mol.parse(ASPIRIN))["sa_score"])


# ------------------------------------------------------------------ efficiency metrics
def test_ligand_efficiency_from_pic50():
    d = props.descriptors(mol.parse(IMATINIB))
    eff = props.efficiency(d, 7.9)
    assert eff["ligand_efficiency"] == pytest.approx(0.29, abs=0.01)
    assert eff["lle"] == pytest.approx(3.31, abs=0.05)
    assert eff["lelp"] > 10  # imatinib buys potency with lipophilicity


def test_nanomolar_potency_converts_to_p_value():
    d = props.descriptors(mol.parse(ASPIRIN))
    assert props.efficiency(d, 25.0, "ic50_nm")["p_activity"] == pytest.approx(7.6, abs=0.01)
    assert props.efficiency(d, 1000.0, "ic50_nm")["p_activity"] == pytest.approx(6.0, abs=0.01)


def test_efficiency_rejects_bad_input():
    d = props.descriptors(mol.parse(ASPIRIN))
    assert props.efficiency(d, None) is None
    assert "error" in props.efficiency(d, 0.0, "ic50_nm")
    assert "error" in props.efficiency(d, 7.0, "nonsense")


# ------------------------------------------------------------------ structural alerts
def test_rhodanine_fires_pains_and_brenk():
    m, _ = mol.standardize(mol.parse(RHODANINE))
    found = alerts.alert_summary(m)
    assert found["alert_count"] >= 2
    names = [a["alert"] for a in found["alerts"]]
    assert any("rhod" in n.lower() for n in names)
    assert all(a["atom_indices"] for a in found["alerts"])


def test_catechol_fires_pains():
    m, _ = mol.standardize(mol.parse("O=C(Nc1ccccc1)c1ccc(O)c(O)c1"))
    assert any("catechol" in a["alert"].lower() for a in alerts.alerts(m))


def test_imatinib_is_clean_on_default_catalogs():
    m, _ = mol.standardize(mol.parse(IMATINIB))
    assert alerts.alert_summary(m)["clean"]


def test_marketed_drugs_still_trip_alerts():
    """The point the alerts-review skill makes: an alert is not a rejection."""
    gef, _ = mol.standardize(mol.parse(GEFITINIB))
    assert alerts.alert_summary(gef)["alert_count"] > 0


def test_unknown_catalog_is_rejected():
    with pytest.raises(ValueError):
        alerts.alerts(mol.parse(ASPIRIN), ["NOT_A_CATALOG"])


# ------------------------------------------------------------------ library I/O
def test_read_example_library():
    rows = mol.read_library(LIB)
    assert len(rows) == 25
    assert rows[0]["id"] == "ASPIRIN"
    assert rows[0]["fields"]["source"] == "reference"


def test_bad_smiles_is_reported_not_dropped():
    ok, failed = mol.load_mols(mol.read_library(LIB))
    assert len(ok) == 24
    assert [f["id"] for f in failed] == ["BROKEN_INPUT"]


def test_smi_and_missing_file(tmp_path):
    p = tmp_path / "frags.smi"
    p.write_text(f"{ASPIRIN} ASA\n{CAFFEINE} CAF\n")
    rows = mol.read_library(str(p))
    assert [r["id"] for r in rows] == ["ASA", "CAF"]
    with pytest.raises(FileNotFoundError):
        mol.read_library(str(tmp_path / "nope.csv"))


def test_csv_without_smiles_column_raises(tmp_path):
    p = tmp_path / "bad.csv"
    p.write_text("name,value\na,1\n")
    with pytest.raises(ValueError, match="no SMILES column"):
        mol.read_library(str(p))


def test_round_trip_csv_and_sdf(tmp_path):
    ok, _ = mol.load_mols(mol.read_library(LIB, limit=5))
    csv_out = mol.write_table(ok, str(tmp_path / "out.csv"))
    sdf_out = mol.write_table(ok, str(tmp_path / "out.sdf"))
    assert csv_out["rows"] == 5 and sdf_out["rows"] == 5
    assert len(mol.read_library(str(tmp_path / "out.sdf"))) == 5


# ------------------------------------------------------------------ screening cascade
def test_dedupe_collapses_the_same_molecule():
    ok, _ = mol.load_mols(mol.read_library(LIB))
    unique, dups = library.dedupe(ok)
    assert len(dups) == 1 and dups[0]["duplicate_of"] == "ASPIRIN"
    assert len(unique) == len(ok) - 1


def test_screen_reports_attrition_by_stage():
    ok, _ = mol.load_mols(mol.read_library(LIB))
    report = library.screen(ok, filters={"mw": [250, 550], "clogp": [None, 5]},
                            require_rules=["Lipinski"])
    assert report["input_rows"] == 24
    assert report["duplicates_removed"] == 1
    assert report["passed"] + report["failed"] == report["screened"]
    stages = [s["stage"] for s in report["stages"]]
    assert any("deduplicate" in s for s in stages)
    assert any("property filters" in s for s in stages)
    for stage in report["stages"]:
        assert stage["in"] - stage["removed"] == stage["out"]


def test_flag_mode_keeps_alerted_compounds_reject_mode_drops_them():
    ok, _ = mol.load_mols(mol.read_library(LIB))
    flagged = library.screen(ok, alert_action="flag")
    rejected = library.screen(ok, alert_action="reject")
    assert flagged["passed"] > rejected["passed"]
    assert flagged["alert_flagged"] > 0
    # The documented lesson: rejecting on alerts discards marketed drugs.
    dropped = {r["id"] for r in rejected["results"] if r["outcome"] == "alert"}
    assert {"GEFITINIB", "ERLOTINIB"} <= dropped


def test_property_filter_bounds_are_inclusive_both_ways():
    ok, _ = mol.load_mols(mol.read_library(LIB))
    low = library.screen(ok, filters={"mw": [400, None]})
    high = library.screen(ok, filters={"mw": [None, 400]})
    assert low["passed"] and high["passed"]
    assert low["passed"] + high["passed"] >= low["screened"]  # boundary compounds counted once each


def test_distribution_summary():
    ok, _ = mol.load_mols(mol.read_library(LIB))
    report = library.screen(ok)
    dist = library.summarize_distribution(report["results"])
    assert dist["mw"]["min"] <= dist["mw"]["median"] <= dist["mw"]["max"]
    assert dist["mw"]["n"] == report["screened"]


# ------------------------------------------------------------------ similarity, clusters, cores
def test_identical_molecule_is_similarity_one():
    a, b = mol.parse(GEFITINIB), mol.parse(GEFITINIB)
    assert similarity.similarity(a, b) == 1.0


def test_analogue_is_more_similar_than_unrelated_compound():
    ok, _ = mol.load_mols(mol.read_library(LIB))
    by_id = {r["id"]: r["mol"] for r in ok}
    close = similarity.similarity(by_id["GEFITINIB"], by_id["GEFITINIB_ANALOGUE"])
    far = similarity.similarity(by_id["GEFITINIB"], by_id["CAFFEINE"])
    assert close > 0.7 > far


def test_search_respects_threshold_and_orders_by_similarity():
    ok, _ = mol.load_mols(mol.read_library(LIB))
    hits = similarity.search(mol.parse(GEFITINIB), ok, threshold=0.4)
    assert hits[0]["id"] == "GEFITINIB" and hits[0]["similarity"] == 1.0
    assert all(h["similarity"] >= 0.4 for h in hits)
    assert hits == sorted(hits, key=lambda h: h["similarity"], reverse=True)


def test_clustering_groups_analogue_pairs():
    ok, _ = mol.load_mols(mol.read_library(LIB))
    result = similarity.cluster(ok, cutoff=0.35)
    pairs = {frozenset(c["member_ids"]) for c in result["clusters"] if c["size"] == 2}
    assert frozenset({"GEFITINIB", "GEFITINIB_ANALOGUE"}) in pairs
    assert frozenset({"ERLOTINIB", "ERLOTINIB_ANALOGUE"}) in pairs
    assert result["cluster_count"] == len(result["clusters"])


def test_looser_cutoff_gives_fewer_clusters():
    ok, _ = mol.load_mols(mol.read_library(LIB))
    tight = similarity.cluster(ok, cutoff=0.2)["cluster_count"]
    loose = similarity.cluster(ok, cutoff=0.6)["cluster_count"]
    assert loose < tight


def test_murcko_scaffold_strips_side_chains():
    assert similarity.murcko_scaffold(mol.parse("CC(=O)Nc1ccc(O)cc1")) == "c1ccccc1"
    assert similarity.scaffold_counts([
        {"id": "a", "mol": mol.parse("CCO")}])["acyclic_molecules"] == 1


def test_mcs_finds_the_quinazoline_core_of_the_egfr_set():
    ok, _ = mol.load_mols(mol.read_library(LIB))
    egfr = [r["mol"] for r in ok if r["id"].startswith(("GEFITINIB", "ERLOTINIB"))]
    core = similarity.mcs(egfr)
    assert core["atoms"] >= 15 and not core["timed_out"]
    assert mol.parse(GEFITINIB).HasSubstructMatch(Chem.MolFromSmarts(core["mcs_smarts"]))


def test_mcs_needs_two_molecules():
    assert "error" in similarity.mcs([mol.parse(ASPIRIN)])


def test_diverse_subset_size_and_passthrough():
    ok, _ = mol.load_mols(mol.read_library(LIB))
    assert len(similarity.diverse_subset(ok, 5)) == 5
    assert len(similarity.diverse_subset(ok, 999)) == len(ok)


# ------------------------------------------------------------------ depiction
def test_grid_and_single_write_images(tmp_path):
    ok, _ = mol.load_mols(mol.read_library(LIB, limit=4))
    png = draw.grid(ok, str(tmp_path / "g.png"), legend_keys=("id",))
    svg = draw.grid(ok, str(tmp_path / "g.svg"), legend_keys=("id",))
    one = draw.single(ok[0]["mol"], str(tmp_path / "one.png"), highlight_atoms=[0, 1])
    assert Path(png["path"]).stat().st_size > 1000
    assert Path(svg["path"]).read_text().lstrip().startswith("<?xml")
    assert Path(one["path"]).stat().st_size > 1000


def test_grid_aligns_to_a_core(tmp_path):
    ok, _ = mol.load_mols(mol.read_library(LIB))
    egfr = [r for r in ok if r["id"].startswith(("GEFITINIB", "ERLOTINIB"))]
    core = similarity.mcs([r["mol"] for r in egfr])
    assert draw.grid(egfr, str(tmp_path / "s.png"), align_to=core["mcs_smarts"])["aligned"]


def test_grid_rejects_empty_input(tmp_path):
    with pytest.raises(ValueError):
        draw.grid([], str(tmp_path / "x.png"))
