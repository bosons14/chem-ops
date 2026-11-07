#!/usr/bin/env python3
"""chem-ops MCP server: compound screening, property profiling, structural alerts, similarity,
clustering, scaffolds and series cores, built on RDKit.

Runs over stdio. Relative paths resolve against CHEMOPS_PROJECT_DIR, else CLAUDE_PROJECT_DIR, else cwd.
Large result sets are written to a file and summarized in the response rather than returned whole.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

try:  # mcp >= 2.0
    from mcp.server.mcpserver import MCPServer as FastMCP
except ImportError:  # mcp 1.x
    from mcp.server.fastmcp import FastMCP

from chemops import __version__, alerts, draw, library, mol, props, similarity

mcp = FastMCP("chem-ops")

#: Rows returned inline per call. The rest go to the output file.
INLINE_LIMIT = 50


def _base() -> Path:
    return Path(os.environ.get("CHEMOPS_PROJECT_DIR") or os.environ.get("CLAUDE_PROJECT_DIR")
                or os.getcwd())


def _path(p: str) -> str:
    q = Path(p).expanduser()
    return str(q if q.is_absolute() else (_base() / q).resolve())


def _clean(rows: list[dict]) -> list[dict]:
    """Drop RDKit objects so the payload serializes."""
    return [{k: v for k, v in row.items() if k != "mol"} for row in rows]


def _as_list(smiles: list[str] | str) -> list[str]:
    if isinstance(smiles, str):
        return [s for s in (part.strip() for part in smiles.replace("\n", " ").split()) if s]
    return list(smiles)


def _load(path: str, smiles_column: str | None = None, id_column: str | None = None,
          limit: int | None = None, standardize: bool = True):
    rows = mol.read_library(_path(path), smiles_column, id_column, limit)
    ok, failed = mol.load_mols(rows, standardize_mols=standardize)
    return ok, failed


def _rows_from(smiles: list[str] | str | None, library_path: str | None,
               smiles_column: str | None = None, id_column: str | None = None,
               limit: int | None = None):
    """Accept either inline SMILES or a library file, the two ways a chemist has compounds to hand."""
    if library_path:
        return _load(library_path, smiles_column, id_column, limit)
    if not smiles:
        raise ValueError("give either smiles or library_path")
    raw = [{"id": f"mol{i}", "smiles": s, "fields": {}} for i, s in enumerate(_as_list(smiles), 1)]
    return mol.load_mols(raw)


# --------------------------------------------------------------------------- status
@mcp.tool()
def chem_tool_status() -> dict:
    """Report the RDKit version, which structural-alert catalogs and fingerprints are available, the
    filterable descriptors, and whether synthetic accessibility scoring is installed.
    Call this first when unsure what this environment can do."""
    import rdkit

    return {
        "version": __version__,
        "rdkit_version": rdkit.__version__,
        "project_dir": str(_base()),
        "sa_score_available": props.sa_available(),
        "alert_catalogs": alerts.catalog_names(),
        "default_alert_catalogs": list(alerts.DEFAULT_CATALOGS),
        "fingerprints": similarity.FINGERPRINTS,
        "filterable_descriptors": sorted(props.PROPERTY_FILTERS),
        "rule_sets": ["Lipinski Ro5", "Veber", "Egan", "Ghose", "Lead-likeness", "Rule of 3"],
        "supported_input_formats": [".csv", ".tsv", ".smi", ".txt", ".sdf", "any of these .gz"],
    }


# --------------------------------------------------------------------------- single compounds
@mcp.tool()
def molecule_profile(smiles: list[str] | str, potency: float | None = None,
                     potency_kind: str = "pic50", include_alerts: bool = True,
                     alert_catalogs: list[str] | None = None,
                     canonical_tautomer: bool = False, image_path: str | None = None) -> dict:
    """Full profile of one or a few compounds: descriptors (MW, cLogP, TPSA, HBD/HBA, rotatable bonds,
    aromatic rings, fraction Csp3, stereocentres, QED, synthetic accessibility), every drug-likeness
    rule set with the specific criteria failed, and structural alerts.

    Supply `potency` to also get ligand efficiency metrics (LE, LLE, LELP): either a p-value
    (potency_kind pic50/pki/pkd) or a nanomolar concentration (ic50_nm/ki_nm/kd_nm).
    `image_path` writes a depiction (.png or .svg) with any alert atoms highlighted."""
    out = []
    for smi in _as_list(smiles):
        m = mol.parse(smi)
        if m is None:
            out.append({"input": smi, "error": "SMILES would not parse"})
            continue
        std, changes = mol.standardize(m, canonical_tautomer=canonical_tautomer)
        item = props.profile(std, potency, potency_kind)
        item["input"] = smi
        item["inchikey"] = mol.inchikey(std)
        if changes:
            item["standardization"] = changes
        if include_alerts:
            item["alerts"] = alerts.alert_summary(std, alert_catalogs)
        out.append(item)

    result = {"molecules": out, "count": len(out)}
    if image_path:
        rows = []
        for item in out:
            if "error" in item:
                continue
            atoms = [i for a in item.get("alerts", {}).get("alerts", []) for i in a["atom_indices"]]
            rows.append({"id": item["input"], "smiles": item["smiles"], "alert_atoms": atoms})
        if rows:
            result["image"] = draw.grid(rows, _path(image_path), legend_keys=("id",),
                                        highlight_key="alert_atoms")
    result["note"] = ("Rule sets are triage guidelines, not verdicts - oral kinase inhibitors, "
                      "macrocycles and antibiotics routinely break them.")
    return result


@mcp.tool()
def structural_alerts(smiles: list[str] | str | None = None, library_path: str | None = None,
                      catalogs: list[str] | None = None, smiles_column: str | None = None,
                      id_column: str | None = None, out_path: str | None = None,
                      image_path: str | None = None, limit: int | None = None) -> dict:
    """Run structural-alert catalogs (PAINS, Brenk, NIH, ZINC, ChEMBL house filters) over compounds
    given inline or in a file, reporting which alert fired and on which atoms.

    Call chem_tool_status for the catalog list. Defaults to PAINS + Brenk + NIH. An alert is a flag
    for a chemist to judge, not a rejection: several marketed drugs carry one."""
    rows, failed = _rows_from(smiles, library_path, smiles_column, id_column, limit)
    catalog_list = catalogs if catalogs is not None else list(alerts.DEFAULT_CATALOGS)

    flagged, clean, histogram = [], 0, {}
    for row in rows:
        info = alerts.alert_summary(row["mol"], catalog_list)
        if info["clean"]:
            clean += 1
            continue
        for item in info["alerts"]:
            histogram[item["alert"]] = histogram.get(item["alert"], 0) + 1
        flagged.append({
            "id": row["id"], "smiles": row["smiles"], "alert_count": info["alert_count"],
            "alerts": [a["alert"] for a in info["alerts"]],
            "alert_atoms": sorted({i for a in info["alerts"] for i in a["atom_indices"]}),
            "mol": row["mol"],
        })

    result = {
        "screened": len(rows),
        "clean": clean,
        "flagged": len(flagged),
        "catalogs": catalog_list,
        "parse_failures": failed[:20],
        "most_common_alerts": sorted(histogram.items(), key=lambda kv: kv[1], reverse=True)[:20],
        "flagged_compounds": _clean(flagged[:INLINE_LIMIT]),
        "truncated": len(flagged) > INLINE_LIMIT,
    }
    if out_path and flagged:
        result["written"] = mol.write_table(flagged, _path(out_path))
    if image_path and flagged:
        result["image"] = draw.grid(flagged[:24], _path(image_path),
                                    legend_keys=("id", "alert_count"), highlight_key="alert_atoms")
    return result


# --------------------------------------------------------------------------- library screening
@mcp.tool()
def screen_library(library_path: str, filters: dict | None = None,
                   require_rules: list[str] | None = None, alert_catalogs: list[str] | None = None,
                   alert_action: str = "flag", max_alerts: int = 0, deduplicate: bool = True,
                   out_passed: str | None = None, out_failed: str | None = None,
                   smiles_column: str | None = None, id_column: str | None = None,
                   limit: int | None = None) -> dict:
    """Screen a compound file through a property/rule/alert cascade and report the attrition at each
    stage, the most common reasons for failure, and the property distribution of what survived.

    filters: {descriptor: [low, high]} with either bound null, e.g.
      {"mw": [250, 450], "clogp": [null, 4], "tpsa": [40, 120], "sa_score": [null, 4.5]}.
      Descriptor names come from chem_tool_status.
    require_rules: rule sets that must pass, e.g. ["Lipinski", "Veber"].
    alert_action: "flag" computes alerts but keeps compounds (recommended - look before discarding);
      "reject" fails anything with more than max_alerts alerts.
    out_passed / out_failed: write the full results to .csv or .sdf (the response is summarized)."""
    rows, parse_failed = _load(library_path, smiles_column, id_column, limit)
    report = library.screen(rows, filters, require_rules, alert_catalogs, alert_action,
                            deduplicate, max_alerts)

    passed = [r for r in report["results"] if r["outcome"] == "pass"]
    failed = [r for r in report["results"] if r["outcome"] != "pass"]
    report["parse_failures"] = parse_failed[:20]
    report["parse_failure_count"] = len(parse_failed)
    report["distribution_passed"] = library.summarize_distribution(passed)
    report["distribution_all"] = library.summarize_distribution(report["results"])

    written = {}
    if out_passed:
        written["passed"] = mol.write_table(passed, _path(out_passed))
    if out_failed:
        written["failed"] = mol.write_table(failed, _path(out_failed))
    if written:
        report["written"] = written

    report["passed_compounds"] = _clean(passed[:INLINE_LIMIT])
    report["failed_examples"] = _clean(failed[:10])
    report["truncated"] = len(passed) > INLINE_LIMIT
    del report["results"]
    if report["pass_rate"] < 0.05 and report["screened"] > 20:
        report["warning"] = ("Under 5% survived. Check the thresholds against the distribution before "
                             "trusting this cut - a near-empty output usually means a mis-set filter.")
    return report


@mcp.tool()
def standardize_library(library_path: str, out_path: str, canonical_tautomer: bool = False,
                        deduplicate: bool = True, smiles_column: str | None = None,
                        id_column: str | None = None, limit: int | None = None) -> dict:
    """Clean a compound file for screening: normalize functional groups, strip salts and solvents to
    the parent, neutralize charges, optionally canonicalize tautomers, and deduplicate by InChIKey.
    Writes a .csv or .sdf and reports what changed, including which compounds would not parse."""
    raw = mol.read_library(_path(library_path), smiles_column, id_column, limit)
    rows, failed = mol.load_mols(raw, standardize_mols=False)

    changed = 0
    for row in rows:
        std, changes = mol.standardize(row["mol"], canonical_tautomer=canonical_tautomer)
        if changes:
            changed += 1
            row["standardization"] = "; ".join(changes)
        row["mol"] = std
        row["original_smiles"] = row["smiles"]
        row["smiles"] = mol.canonical(std)
        row["inchikey"] = mol.inchikey(std)

    duplicates: list[dict] = []
    if deduplicate:
        rows, duplicates = library.dedupe(rows)

    return {
        "input_rows": len(raw),
        "parsed": len(rows) + len(duplicates),
        "parse_failures": failed[:20],
        "parse_failure_count": len(failed),
        "standardized": changed,
        "duplicates_removed": len(duplicates),
        "duplicates": duplicates[:50],
        "output_rows": len(rows),
        "written": mol.write_table(rows, _path(out_path)),
    }


# --------------------------------------------------------------------------- series and neighbours
@mcp.tool()
def similarity_search(query_smiles: str, library_path: str, threshold: float = 0.7,
                      fingerprint: str = "morgan2", top_n: int | None = 100,
                      smiles_column: str | None = None, id_column: str | None = None,
                      out_path: str | None = None, image_path: str | None = None,
                      limit: int | None = None) -> dict:
    """Find compounds in a file similar to a query structure, by Tanimoto similarity on fingerprints.

    Use this to pull the neighbours of a hit out of a deck before ordering analogues. Treat the score
    as relative: 0.85 on ECFP4 is the common "same chemotype" cut, but it depends on the series."""
    query = mol.parse(query_smiles)
    if query is None:
        return {"error": f"query SMILES would not parse: {query_smiles!r}"}
    query, _ = mol.standardize(query)
    rows, failed = _load(library_path, smiles_column, id_column, limit)
    hits = similarity.search(query, rows, threshold, fingerprint, top_n)

    by_id = {r["id"]: r for r in rows}
    hit_rows = [dict(h, mol=by_id[h["id"]]["mol"]) for h in hits if h["id"] in by_id]
    result = {
        "query": mol.canonical(query),
        "searched": len(rows),
        "hits": len(hits),
        "threshold": threshold,
        "fingerprint": fingerprint,
        "parse_failure_count": len(failed),
        "nearest": hits[:INLINE_LIMIT],
        "truncated": len(hits) > INLINE_LIMIT,
    }
    if out_path and hit_rows:
        result["written"] = mol.write_table(hit_rows, _path(out_path))
    if image_path and hit_rows:
        result["image"] = draw.grid(hit_rows[:24], _path(image_path),
                                    legend_keys=("id", "similarity"))
    return result


@mcp.tool()
def cluster_library(library_path: str, cutoff: float = 0.35, fingerprint: str = "morgan2",
                    pick_per_cluster: int = 1, smiles_column: str | None = None,
                    id_column: str | None = None, out_path: str | None = None,
                    image_path: str | None = None, limit: int | None = None) -> dict:
    """Butina-cluster a compound set into chemotypes and return the clusters largest first with their
    centroids - the representative set to test or buy.

    cutoff is a Tanimoto *distance*: 0.35 means 0.65 similarity, the usual chemotype grouping on
    ECFP4. Lower cutoff gives tighter, more numerous clusters. out_path writes one row per compound
    with its cluster number; image_path draws the cluster centroids."""
    rows, failed = _load(library_path, smiles_column, id_column, limit)
    if not rows:
        return {"error": "no parsable compounds in file", "parse_failures": failed[:20]}
    result = similarity.cluster(rows, cutoff, fingerprint)

    by_id = {r["id"]: r for r in rows}
    assignment = []
    for cluster_info in result["clusters"]:
        for rank, member_id in enumerate(cluster_info["member_ids"]):
            row = by_id.get(member_id)
            if row is None:
                continue
            assignment.append({
                "id": member_id, "smiles": row["smiles"], "cluster": cluster_info["cluster"],
                "cluster_size": cluster_info["size"], "is_centroid": rank == 0, "mol": row["mol"],
            })

    representatives = []
    for cluster_info in result["clusters"]:
        members = [by_id[m] for m in cluster_info["member_ids"] if m in by_id]
        chosen = members[:1] if pick_per_cluster <= 1 else similarity.diverse_subset(
            members, min(pick_per_cluster, len(members)), fingerprint)
        for row in chosen:
            representatives.append({"id": row["id"], "smiles": row["smiles"],
                                    "cluster": cluster_info["cluster"],
                                    "cluster_size": cluster_info["size"], "mol": row["mol"]})

    out = {
        "clustered": len(rows),
        "parse_failure_count": len(failed),
        "cluster_count": result["cluster_count"],
        "singletons": result["singletons"],
        "cutoff_distance": cutoff,
        "cutoff_similarity": result["cutoff_similarity"],
        "fingerprint": fingerprint,
        "largest_clusters": [
            {k: v for k, v in c.items() if k != "member_ids"} | {"member_ids": c["member_ids"][:20]}
            for c in result["clusters"][:INLINE_LIMIT]
        ],
        "representatives": _clean(representatives[:INLINE_LIMIT]),
        "representative_count": len(representatives),
    }
    if out_path and assignment:
        out["written"] = mol.write_table(assignment, _path(out_path))
    if image_path and representatives:
        out["image"] = draw.grid(representatives[:24], _path(image_path),
                                 legend_keys=("id", "cluster", "cluster_size"))
    return out


@mcp.tool()
def scaffold_analysis(library_path: str | None = None, smiles: list[str] | str | None = None,
                      generic: bool = False, top_n: int = 25, smiles_column: str | None = None,
                      id_column: str | None = None, out_path: str | None = None,
                      image_path: str | None = None, limit: int | None = None) -> dict:
    """Group compounds by Bemis-Murcko scaffold and report the most frequent scaffolds with their
    members - how many distinct chemotypes a hit list really contains.

    generic=true reduces scaffolds to the ring skeleton (atom and bond types erased), which merges
    bioisosteric variants of one core."""
    rows, failed = _rows_from(smiles, library_path, smiles_column, id_column, limit)
    result = similarity.scaffold_counts(rows, generic, top_n)
    result["analyzed"] = len(rows)
    result["parse_failure_count"] = len(failed)

    scaffold_rows = [{"id": f"scaffold{i}", "smiles": s["scaffold"], "count": s["count"]}
                     for i, s in enumerate(result["scaffolds"], 1) if s["scaffold"]]
    if out_path:
        by_id = {r["id"]: r for r in rows}
        assignment = []
        for item in result["scaffolds"]:
            for member in item["member_ids"]:
                row = by_id.get(member)
                if row:
                    assignment.append({"id": member, "smiles": row["smiles"],
                                       "scaffold": item["scaffold"], "scaffold_count": item["count"]})
        result["written"] = mol.write_table(assignment, _path(out_path))
    if image_path and scaffold_rows:
        result["image"] = draw.grid(scaffold_rows[:24], _path(image_path),
                                    legend_keys=("id", "count"))
    return result


@mcp.tool()
def series_core(smiles: list[str] | str | None = None, library_path: str | None = None,
                threshold: float = 1.0, timeout: int = 20, smiles_column: str | None = None,
                id_column: str | None = None, image_path: str | None = None,
                limit: int | None = None) -> dict:
    """Find the maximum common substructure of a set of compounds - the shared core of a series, and
    the starting point for an SAR table or an R-group decomposition.

    threshold below 1.0 finds a core present in that fraction of the set, which is what you want when
    one or two members are structurally off. image_path draws the series aligned on the core."""
    rows, failed = _rows_from(smiles, library_path, smiles_column, id_column, limit)
    if len(rows) < 2:
        return {"error": "need at least 2 parsable compounds", "parse_failures": failed[:20]}
    result = similarity.mcs([r["mol"] for r in rows], timeout=timeout, threshold=threshold)
    result["parse_failure_count"] = len(failed)
    if image_path and result.get("mcs_smarts"):
        result["image"] = draw.grid(rows[:24], _path(image_path), legend_keys=("id",),
                                    align_to=result["mcs_smarts"])
    return result


@mcp.tool()
def depict_molecules(smiles: list[str] | str | None = None, library_path: str | None = None,
                     image_path: str = "molecules.png", per_row: int = 4, size: int = 300,
                     align_to: str | None = None, legend_fields: list[str] | None = None,
                     smiles_column: str | None = None, id_column: str | None = None,
                     limit: int | None = 24) -> dict:
    """Draw compounds as a 2D grid image (.png or .svg) so they can actually be looked at.

    align_to takes a SMARTS or SMILES core and orients every structure on it, which is how a series
    should be presented. legend_fields names columns from the input file to print under each structure."""
    rows, failed = _rows_from(smiles, library_path, smiles_column, id_column, limit)
    if not rows:
        return {"error": "no parsable compounds", "parse_failures": failed[:20]}
    for row in rows:
        for key, value in (row.get("fields") or {}).items():
            row.setdefault(key, value)
    keys = tuple(legend_fields) if legend_fields else ("id",)
    result = draw.grid(rows, _path(image_path), per_row=per_row, size=size,
                       legend_keys=keys, align_to=align_to)
    result["parse_failure_count"] = len(failed)
    return result


if __name__ == "__main__":
    mcp.run()
