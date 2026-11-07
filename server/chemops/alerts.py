"""Structural alerts from RDKit's filter catalogs (PAINS, Brenk, NIH, ZINC, ChEMBL house filters).

An alert is a *flag for a chemist to judge*, not a rejection. PAINS in particular were derived from
AlphaScreen output: a PAINS hit in a biochemical assay deserves orthogonal confirmation, while the
same substructure in a series with clean SAR and a crystal structure may be entirely fine. Several
marketed drugs carry Brenk or PAINS alerts.
"""

from __future__ import annotations

from functools import lru_cache

from rdkit import Chem
from rdkit.Chem import FilterCatalog
from rdkit.Chem.FilterCatalog import FilterCatalogParams

#: Catalog name -> (RDKit enum, what it is for).
CATALOGS = {
    "PAINS": (FilterCatalogParams.FilterCatalogs.PAINS,
              ("Pan-assay interference compounds (Baell & Holloway) - frequent hitters in "
               "biochemical screens")),
    "PAINS_A": (FilterCatalogParams.FilterCatalogs.PAINS_A, "PAINS family A - strongest evidence"),
    "PAINS_B": (FilterCatalogParams.FilterCatalogs.PAINS_B, "PAINS family B"),
    "PAINS_C": (FilterCatalogParams.FilterCatalogs.PAINS_C, "PAINS family C - weakest evidence"),
    "BRENK": (FilterCatalogParams.FilterCatalogs.BRENK,
              "Brenk unwanted functionality - reactive, metabolically labile or toxicophoric groups"),
    "NIH": (FilterCatalogParams.FilterCatalogs.NIH,
            "NIH filters - reactive and assay-interfering functionality"),
    "ZINC": (FilterCatalogParams.FilterCatalogs.ZINC,
             "ZINC drug-likeness and reactivity filters"),
    "CHEMBL_BMS": (FilterCatalogParams.FilterCatalogs.CHEMBL_BMS, "BMS HTS deck filters"),
    "CHEMBL_Dundee": (FilterCatalogParams.FilterCatalogs.CHEMBL_Dundee, "University of Dundee filters"),
    "CHEMBL_Glaxo": (FilterCatalogParams.FilterCatalogs.CHEMBL_Glaxo, "Glaxo hard filters"),
    "CHEMBL_Inpharmatica": (FilterCatalogParams.FilterCatalogs.CHEMBL_Inpharmatica,
                            "Inpharmatica unwanted fragments"),
    "CHEMBL_LINT": (FilterCatalogParams.FilterCatalogs.CHEMBL_LINT, "Pfizer LINT filters"),
    "CHEMBL_MLSMR": (FilterCatalogParams.FilterCatalogs.CHEMBL_MLSMR, "NIH MLSMR excluded functionality"),
    "CHEMBL_SureChEMBL": (FilterCatalogParams.FilterCatalogs.CHEMBL_SureChEMBL,
                          "SureChEMBL non-MedChem-friendly filters"),
}

DEFAULT_CATALOGS = ("PAINS", "BRENK", "NIH")


def catalog_names() -> dict:
    return {name: meaning for name, (_, meaning) in CATALOGS.items()}


@lru_cache(maxsize=16)
def _catalog(names: tuple[str, ...]) -> FilterCatalog.FilterCatalog:
    params = FilterCatalogParams()
    for name in names:
        params.AddCatalog(CATALOGS[name][0])
    return FilterCatalog.FilterCatalog(params)


def _resolve(catalogs: list[str] | tuple[str, ...] | None) -> tuple[str, ...]:
    names = tuple(catalogs) if catalogs else DEFAULT_CATALOGS
    unknown = [n for n in names if n not in CATALOGS]
    if unknown:
        raise ValueError(f"unknown catalog(s) {unknown}; available: {sorted(CATALOGS)}")
    return names


def alerts(mol: Chem.Mol, catalogs: list[str] | None = None) -> list[dict]:
    """Every alert firing on this molecule, with the matching atoms so a chemist can see where."""
    names = _resolve(catalogs)
    catalog = _catalog(names)
    found = []
    for match in catalog.GetMatches(mol):
        available = set(match.GetPropList())
        source = match.GetProp("FilterSet") if "FilterSet" in available else ""
        scope = match.GetProp("Scope") if "Scope" in available else ""
        atoms: list[int] = []
        for filter_match in match.GetFilterMatches(mol):
            atoms.extend(int(pair[1]) for pair in filter_match.atomPairs)
        found.append({
            "alert": match.GetDescription(),
            "catalog": source or ",".join(names),
            "why": scope,
            "atom_indices": sorted(set(atoms)),
        })
    # One substructure can be listed in several catalogs; keep each distinct alert once.
    seen, unique = set(), []
    for item in found:
        key = (item["alert"], tuple(item["atom_indices"]))
        if key not in seen:
            seen.add(key)
            unique.append(item)
    return unique


def alert_summary(mol: Chem.Mol, catalogs: list[str] | None = None) -> dict:
    found = alerts(mol, catalogs)
    by_catalog: dict[str, int] = {}
    for item in found:
        by_catalog[item["catalog"]] = by_catalog.get(item["catalog"], 0) + 1
    return {
        "alert_count": len(found),
        "alerts": found,
        "by_catalog": by_catalog,
        "clean": not found,
    }
