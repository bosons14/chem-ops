"""Physicochemical descriptors, drug-likeness rules and efficiency metrics.

The rule sets here are the published filters medicinal chemists use to triage screening output. They
are *guidelines*, not pass/fail verdicts: every marketed drug class has members that break them
(macrocycles, antibiotics, kinase inhibitors beyond Ro5). Report them as such.
"""

from __future__ import annotations

import math
import os
import sys

from rdkit import Chem
from rdkit.Chem import QED, Crippen, Descriptors, Lipinski, RDConfig, rdMolDescriptors

# Synthetic accessibility lives in RDKit's Contrib tree, which is not importable by default.
_SA_AVAILABLE = False
try:
    sys.path.append(os.path.join(RDConfig.RDContribDir, "SA_Score"))
    import sascorer

    _SA_AVAILABLE = True
except Exception:  # noqa: BLE001 # pragma: no cover - depends on the RDKit build
    sascorer = None


def sa_available() -> bool:
    return _SA_AVAILABLE


# --------------------------------------------------------------------------- descriptors
def descriptors(mol: Chem.Mol) -> dict:
    """The descriptor set that drives screening decisions, rounded to the precision worth reporting."""
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    heavy = mol.GetNumHeavyAtoms()
    carbons = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 6)
    d = {
        "formula": rdMolDescriptors.CalcMolFormula(mol),
        "mw": round(Descriptors.MolWt(mol), 2),
        "heavy_atoms": heavy,
        "clogp": round(Crippen.MolLogP(mol), 2),
        "tpsa": round(rdMolDescriptors.CalcTPSA(mol), 2),
        "hbd": Lipinski.NumHDonors(mol),
        "hba": Lipinski.NumHAcceptors(mol),
        "rotatable_bonds": Lipinski.NumRotatableBonds(mol),
        "rings": ring_info.NumRings(),
        "aromatic_rings": rdMolDescriptors.CalcNumAromaticRings(mol),
        "largest_ring": max((len(r) for r in rings), default=0),
        "fraction_csp3": round(rdMolDescriptors.CalcFractionCSP3(mol), 3),
        "heteroatoms": rdMolDescriptors.CalcNumHeteroatoms(mol),
        "formal_charge": Chem.GetFormalCharge(mol),
        "stereocenters": rdMolDescriptors.CalcNumAtomStereoCenters(mol),
        "unspecified_stereocenters": rdMolDescriptors.CalcNumUnspecifiedAtomStereoCenters(mol),
        "molar_refractivity": round(Crippen.MolMR(mol), 2),
        "qed": round(QED.qed(mol), 3),
        "heavy_atom_logp_ratio": round(Crippen.MolLogP(mol) / heavy, 3) if heavy else None,
        "carbon_fraction": round(carbons / heavy, 3) if heavy else None,
    }
    if _SA_AVAILABLE:
        d["sa_score"] = round(sascorer.calculateScore(mol), 2)
    return d


# --------------------------------------------------------------------------- rule sets
def _rule(name: str, checks: list[tuple[str, bool]], allowed_violations: int = 0) -> dict:
    failed = [label for label, ok in checks if not ok]
    return {
        "rule": name,
        "violations": len(failed),
        "failed_criteria": failed,
        "pass": len(failed) <= allowed_violations,
    }


def rule_sets(mol: Chem.Mol, d: dict | None = None) -> list[dict]:
    """Evaluate the standard filters. ``pass`` follows each rule's own published tolerance: Lipinski
    allows one violation, the others allow none."""
    d = d or descriptors(mol)
    mw, logp, tpsa = d["mw"], d["clogp"], d["tpsa"]
    hbd, hba, rotb = d["hbd"], d["hba"], d["rotatable_bonds"]
    heavy, mr = d["heavy_atoms"], d["molar_refractivity"]

    out = [
        _rule("Lipinski Ro5", [
            (f"MW {mw} > 500", mw <= 500),
            (f"cLogP {logp} > 5", logp <= 5),
            (f"HBD {hbd} > 5", hbd <= 5),
            (f"HBA {hba} > 10", hba <= 10),
        ], allowed_violations=1),
        _rule("Veber (oral bioavailability)", [
            (f"rotatable bonds {rotb} > 10", rotb <= 10),
            (f"TPSA {tpsa} > 140", tpsa <= 140),
        ]),
        _rule("Egan (absorption)", [
            (f"TPSA {tpsa} > 131.6", tpsa <= 131.6),
            (f"cLogP {logp} > 5.88", logp <= 5.88),
        ]),
        _rule("Ghose", [
            (f"MW {mw} outside 160-480", 160 <= mw <= 480),
            (f"cLogP {logp} outside -0.4-5.6", -0.4 <= logp <= 5.6),
            (f"heavy atoms {heavy} outside 20-70", 20 <= heavy <= 70),
            (f"molar refractivity {mr} outside 40-130", 40 <= mr <= 130),
        ]),
        _rule("Lead-likeness", [
            (f"MW {mw} > 350", mw <= 350),
            (f"cLogP {logp} > 3.5", logp <= 3.5),
            (f"rotatable bonds {rotb} > 7", rotb <= 7),
        ]),
        _rule("Rule of 3 (fragments)", [
            (f"MW {mw} > 300", mw <= 300),
            (f"cLogP {logp} > 3", logp <= 3),
            (f"HBD {hbd} > 3", hbd <= 3),
            (f"HBA {hba} > 3", hba <= 3),
            (f"rotatable bonds {rotb} > 3", rotb <= 3),
        ]),
    ]
    return out


# --------------------------------------------------------------------------- efficiency metrics
def efficiency(d: dict, potency: float | None = None,
               potency_kind: str = "pic50") -> dict | None:
    """Ligand efficiency metrics, which need measured potency and so are only computed when given.

    ``potency_kind``: ``pic50``/``pki``/``pkd`` (already -log10 M), or ``ic50_nm``/``ki_nm``/``kd_nm``
    in nanomolar. LE is kcal/mol per heavy atom; LLE is pActivity - cLogP; LELP is cLogP / LE.
    """
    if potency is None:
        return None
    kind = potency_kind.lower()
    if kind in ("ic50_nm", "ki_nm", "kd_nm"):
        if potency <= 0:
            return {"error": "potency in nM must be greater than 0"}
        p = -math.log10(potency * 1e-9)
    elif kind in ("pic50", "pki", "pkd"):
        p = potency
    else:
        return {"error": f"unknown potency_kind {potency_kind!r}"}

    heavy = d.get("heavy_atoms") or 0
    if not heavy:
        return {"error": "no heavy atoms"}
    le = 1.37 * p / heavy  # 1.37 = 2.303 RT at 298 K, kcal/mol
    lle = p - d["clogp"]
    return {
        "p_activity": round(p, 2),
        "ligand_efficiency": round(le, 3),
        "lle": round(lle, 2),
        "lelp": round(d["clogp"] / le, 1) if le else None,
        "guide": "LE >= 0.3 and LLE >= 5 are the usual optimization targets; LELP above 10 flags "
                 "potency bought with lipophilicity.",
    }


# --------------------------------------------------------------------------- profile
def profile(mol: Chem.Mol, potency: float | None = None, potency_kind: str = "pic50") -> dict:
    """Descriptors plus every rule set, and efficiency metrics when potency is supplied."""
    d = descriptors(mol)
    rules = rule_sets(mol, d)
    out = {
        "smiles": Chem.MolToSmiles(mol),
        "descriptors": d,
        "rules": rules,
        "rules_failed": [r["rule"] for r in rules if not r["pass"]],
    }
    eff = efficiency(d, potency, potency_kind)
    if eff:
        out["efficiency"] = eff
    return out


# --------------------------------------------------------------------------- filter predicates
#: Named numeric filters usable in a screening cascade, as (descriptor key, low, high).
PROPERTY_FILTERS = {
    "mw": ("mw", 0.0, 500.0),
    "clogp": ("clogp", -float("inf"), 5.0),
    "tpsa": ("tpsa", 0.0, 140.0),
    "hbd": ("hbd", 0, 5),
    "hba": ("hba", 0, 10),
    "rotatable_bonds": ("rotatable_bonds", 0, 10),
    "heavy_atoms": ("heavy_atoms", 0, 70),
    "qed": ("qed", 0.0, 1.0),
    "sa_score": ("sa_score", 1.0, 6.0),
    "aromatic_rings": ("aromatic_rings", 0, 5),
    "fraction_csp3": ("fraction_csp3", 0.0, 1.0),
    "formal_charge": ("formal_charge", -2, 2),
}
