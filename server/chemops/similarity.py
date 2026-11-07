"""Fingerprints, similarity search, clustering, scaffolds and maximum common substructure.

These are the operations that turn a flat hit list into series a chemist can reason about: which hits
are the same chemotype, what the common core is, and which compound best represents each cluster.
"""

from __future__ import annotations

from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator, rdFMCS
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.ML.Cluster import Butina

FINGERPRINTS = {
    "morgan2": "ECFP4-equivalent Morgan, radius 2, 2048 bits - the default for similarity and clustering",
    "morgan3": "ECFP6-equivalent Morgan, radius 3",
    "rdkit": "RDKit path-based fingerprint",
    "atompair": "Atom-pair fingerprint - better for small or very diverse sets",
    "maccs": "166 MACCS structural keys - coarse, fast, good for a first pass",
}


def _generator(kind: str):
    kind = kind.lower()
    if kind in ("morgan2", "ecfp4", "morgan"):
        return rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    if kind in ("morgan3", "ecfp6"):
        return rdFingerprintGenerator.GetMorganGenerator(radius=3, fpSize=2048)
    if kind == "rdkit":
        return rdFingerprintGenerator.GetRDKitFPGenerator(fpSize=2048)
    if kind in ("atompair", "ap"):
        return rdFingerprintGenerator.GetAtomPairGenerator(fpSize=2048)
    raise ValueError(f"unknown fingerprint {kind!r}; available: {sorted(FINGERPRINTS)}")


def fingerprints(mols: list[Chem.Mol], kind: str = "morgan2") -> list:
    if kind.lower() == "maccs":
        from rdkit.Chem import MACCSkeys
        return [MACCSkeys.GenMACCSKeys(m) for m in mols]
    gen = _generator(kind)
    return [gen.GetFingerprint(m) for m in mols]


def similarity(a: Chem.Mol, b: Chem.Mol, kind: str = "morgan2") -> float:
    fa, fb = fingerprints([a, b], kind)
    return round(DataStructs.TanimotoSimilarity(fa, fb), 4)


def search(query: Chem.Mol, rows: list[dict], threshold: float = 0.7,
           kind: str = "morgan2", top_n: int | None = None) -> list[dict]:
    """Tanimoto similarity of every row against the query, above ``threshold``, best first.

    Tanimoto on ECFP4 is not a distance a chemist should over-read: 0.85 is often quoted as "similar"
    but the right cut depends on the fingerprint and the series. Returned scores let the user judge.
    """
    mols = [r["mol"] for r in rows]
    fps = fingerprints([query] + mols, kind)
    qfp, rest = fps[0], fps[1:]
    scores = DataStructs.BulkTanimotoSimilarity(qfp, rest)
    hits = [
        {"id": row["id"], "smiles": row["smiles"], "similarity": round(score, 4)}
        for row, score in zip(rows, scores) if score >= threshold
    ]
    hits.sort(key=lambda h: h["similarity"], reverse=True)
    return hits[:top_n] if top_n else hits


def cluster(rows: list[dict], cutoff: float = 0.35, kind: str = "morgan2") -> dict:
    """Butina clustering on Tanimoto distance. ``cutoff`` is a *distance*: 0.35 distance is 0.65
    similarity, the usual chemotype-level grouping for ECFP4.

    Returns clusters ordered largest first, each with its centroid (the Butina cluster head).
    """
    if not rows:
        return {"clusters": [], "cluster_count": 0, "singletons": 0}
    fps = fingerprints([r["mol"] for r in rows], kind)
    n = len(fps)
    dists: list[float] = []
    for i in range(1, n):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i])
        dists.extend(1 - s for s in sims)
    raw = Butina.ClusterData(dists, n, cutoff, isDistData=True)

    clusters = []
    for idx, members in enumerate(sorted(raw, key=len, reverse=True), start=1):
        centroid = rows[members[0]]
        clusters.append({
            "cluster": idx,
            "size": len(members),
            "centroid_id": centroid["id"],
            "centroid_smiles": centroid["smiles"],
            "member_ids": [rows[m]["id"] for m in members],
        })
    return {
        "clusters": clusters,
        "cluster_count": len(clusters),
        "singletons": sum(1 for c in clusters if c["size"] == 1),
        "cutoff_distance": cutoff,
        "cutoff_similarity": round(1 - cutoff, 3),
        "fingerprint": kind,
    }


def murcko_scaffold(mol: Chem.Mol, generic: bool = False) -> str:
    """Bemis-Murcko scaffold SMILES. ``generic`` strips atom and bond types to the ring skeleton."""
    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    if generic:
        scaffold = MurckoScaffold.MakeScaffoldGeneric(scaffold)
    return Chem.MolToSmiles(scaffold)


def scaffold_counts(rows: list[dict], generic: bool = False, top_n: int = 25) -> dict:
    """Frequency of each Bemis-Murcko scaffold, most common first, with the members of each."""
    groups: dict[str, list[str]] = {}
    acyclic = 0
    for row in rows:
        try:
            key = murcko_scaffold(row["mol"], generic)
        except Exception:  # noqa: BLE001,S112 - one unscaffoldable molecule must not stop the set
            continue
        if not key:  # acyclic molecule: no scaffold
            acyclic += 1
            continue
        groups.setdefault(key, []).append(row["id"])
    ordered = sorted(groups.items(), key=lambda kv: len(kv[1]), reverse=True)
    return {
        "scaffold_count": len(ordered),
        "acyclic_molecules": acyclic,
        "generic": generic,
        "scaffolds": [
            {"scaffold": smi, "count": len(ids), "member_ids": ids[:50]}
            for smi, ids in ordered[:top_n]
        ],
    }


def mcs(mols: list[Chem.Mol], timeout: int = 20, ring_matches_ring_only: bool = True,
        complete_rings_only: bool = True, threshold: float = 1.0) -> dict:
    """Maximum common substructure of a set - the common core of a series.

    ``threshold`` below 1.0 finds a core present in that fraction of the molecules, which is what you
    want when one or two members of a series are structurally off.
    """
    if len(mols) < 2:
        return {"error": "need at least 2 molecules for an MCS"}
    result = rdFMCS.FindMCS(
        mols,
        timeout=timeout,
        ringMatchesRingOnly=ring_matches_ring_only,
        completeRingsOnly=complete_rings_only,
        threshold=threshold,
        atomCompare=rdFMCS.AtomCompare.CompareElements,
        bondCompare=rdFMCS.BondCompare.CompareOrderExact,
    )
    if not result.smartsString:
        return {"error": "no common substructure found", "timed_out": result.canceled}
    core = Chem.MolFromSmarts(result.smartsString)
    smallest = min(m.GetNumHeavyAtoms() for m in mols)
    return {
        "mcs_smarts": result.smartsString,
        "mcs_smiles": Chem.MolToSmiles(core) if core else "",
        "atoms": result.numAtoms,
        "bonds": result.numBonds,
        "fraction_of_smallest": round(result.numAtoms / smallest, 3) if smallest else None,
        "timed_out": result.canceled,
        "molecules": len(mols),
    }


def diverse_subset(rows: list[dict], n: int, kind: str = "morgan2", seed: int = 0xF00D) -> list[dict]:
    """MaxMin diverse pick - the right way to choose ``n`` compounds to buy or test from a cluster set."""
    from rdkit.SimDivFilters import rdSimDivPickers

    if n >= len(rows):
        return rows
    fps = fingerprints([r["mol"] for r in rows], kind)
    picker = rdSimDivPickers.MaxMinPicker()
    idx = picker.LazyBitVectorPick(fps, len(fps), n, seed=seed)
    return [rows[i] for i in idx]
