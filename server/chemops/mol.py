"""Molecule parsing, standardization and library file I/O.

Every public function takes and returns plain data (str/dict/list) so the MCP layer stays thin.
A molecule that will not parse is never dropped silently: it comes back with an ``error`` field.
"""

from __future__ import annotations

import csv
import gzip
import io
import os
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem import inchi
from rdkit.Chem.MolStandardize import rdMolStandardize

# RDKit writes parse errors to stderr, which is the MCP log channel. Collect them per call instead.
RDLogger.DisableLog("rdApp.*")

SMILES_COLUMNS = ("smiles", "smile", "canonical_smiles", "structure", "mol", "smi")
ID_COLUMNS = ("id", "name", "compound_id", "compound", "molecule_id", "mol_id",
              "chembl_id", "catalog_id", "cid", "idnumber", "title")


# --------------------------------------------------------------------------- parsing
def parse(smiles: str, sanitize: bool = True) -> Chem.Mol | None:
    """Parse one SMILES string. Returns None if it will not parse."""
    if not smiles or not smiles.strip():
        return None
    return Chem.MolFromSmiles(smiles.strip(), sanitize=sanitize)


def canonical(mol: Chem.Mol) -> str:
    return Chem.MolToSmiles(mol)


def inchikey(mol: Chem.Mol) -> str:
    """InChIKey, the right identity key for deduplication. Empty string if InChI generation fails."""
    try:
        return inchi.MolToInchiKey(mol) or ""
    except Exception:  # noqa: BLE001 - InChI generation fails in many unrelated ways; no key is fine
        return ""


def standardize(mol: Chem.Mol, strip_salts: bool = True, neutralize: bool = True,
                canonical_tautomer: bool = False) -> tuple[Chem.Mol, list[str]]:
    """Clean a molecule for screening: normalize functional groups, keep the largest organic fragment,
    neutralize charges. Returns the new molecule and the list of changes applied.

    ``canonical_tautomer`` is off by default: it is correct but slow (tens of ms per molecule), so it
    is worth enabling for a single compound profile and not for a 100k library.
    """
    changes: list[str] = []
    before = Chem.MolToSmiles(mol)

    mol = rdMolStandardize.Normalize(mol)
    if Chem.MolToSmiles(mol) != before:
        changes.append("normalized functional groups")

    if strip_salts and len(Chem.GetMolFrags(mol)) > 1:
        mol = rdMolStandardize.FragmentParent(mol)
        changes.append("kept largest fragment (salt/solvent removed)")

    if neutralize:
        charged = Chem.MolToSmiles(mol)
        mol = rdMolStandardize.Uncharger().uncharge(mol)
        if Chem.MolToSmiles(mol) != charged:
            changes.append("neutralized charges")

    if canonical_tautomer:
        taut = Chem.MolToSmiles(mol)
        mol = rdMolStandardize.TautomerEnumerator().Canonicalize(mol)
        if Chem.MolToSmiles(mol) != taut:
            changes.append("canonical tautomer")

    # The standardizer returns molecules whose ring info is not initialized, which makes any later
    # fingerprint or ring descriptor raise. Re-sanitize once here so every caller gets a usable mol.
    try:
        Chem.SanitizeMol(mol)
    except Exception:  # noqa: BLE001 - fall back to the minimum a caller needs, whatever failed
        mol.UpdatePropertyCache(strict=False)
        Chem.FastFindRings(mol)

    return mol, changes


# --------------------------------------------------------------------------- library I/O
def _open_text(path: str):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, encoding="utf-8", errors="replace")


def _pick(fieldnames: list[str], candidates: tuple[str, ...]) -> str | None:
    lowered = {f.lower().strip(): f for f in fieldnames if f}
    for c in candidates:
        if c in lowered:
            return lowered[c]
    for low, orig in lowered.items():  # substring match: "Smiles (canonical)"
        if any(c in low for c in candidates):
            return orig
    return None


def _read_csv(path: str, smiles_col: str | None, id_col: str | None) -> list[dict]:
    with _open_text(path) as fh:
        sample = fh.read(8192)
        fh.seek(0)
        try:
            dialect = csv.Sniffer().sniff(sample, delimiters=",\t;|")
        except csv.Error:
            dialect = csv.excel
        reader = csv.DictReader(fh, dialect=dialect)
        fields = reader.fieldnames or []
        scol = smiles_col or _pick(fields, SMILES_COLUMNS)
        icol = id_col or _pick(fields, ID_COLUMNS)
        if not scol:
            raise ValueError(
                f"no SMILES column found in {Path(path).name}; columns are {fields}. "
                "Pass smiles_column explicitly.")
        rows = []
        for i, row in enumerate(reader, start=1):
            rows.append({
                "id": (row.get(icol) or f"row{i}").strip() if icol else f"row{i}",
                "smiles": (row.get(scol) or "").strip(),
                "fields": {k: v for k, v in row.items() if k not in (scol, icol) and k},
            })
        return rows


def _read_smi(path: str) -> list[dict]:
    rows = []
    with _open_text(path) as fh:
        for i, line in enumerate(fh, start=1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            smi = parts[0]
            if i == 1 and smi.lower() in SMILES_COLUMNS:  # header line
                continue
            rows.append({"id": parts[1] if len(parts) > 1 else f"line{i}", "smiles": smi, "fields": {}})
    return rows


def _read_sdf(path: str) -> list[dict]:
    if str(path).endswith(".gz"):
        with gzip.open(path, "rb") as fh:
            supplier = Chem.ForwardSDMolSupplier(io.BytesIO(fh.read()))
    else:
        supplier = Chem.ForwardSDMolSupplier(path)
    rows = []
    for i, mol in enumerate(supplier, start=1):
        if mol is None:
            rows.append({"id": f"record{i}", "smiles": "", "fields": {},
                         "error": "SDF record would not parse"})
            continue
        props = mol.GetPropsAsDict()
        name = mol.GetProp("_Name") if mol.HasProp("_Name") else ""
        ident = name.strip() or str(props.get("ID") or props.get("id") or f"record{i}")
        rows.append({"id": ident, "smiles": Chem.MolToSmiles(mol),
                     "fields": {k: v for k, v in props.items() if not k.startswith("_")}})
    return rows


def read_library(path: str, smiles_column: str | None = None, id_column: str | None = None,
                 limit: int | None = None) -> list[dict]:
    """Read a compound file into rows of ``{id, smiles, fields}``.

    Supports .csv/.tsv (SMILES and ID columns auto-detected), .smi/.txt (SMILES [ID] per line) and
    .sdf, each optionally .gz. Rows that will not parse keep an ``error`` field rather than vanishing.
    """
    if not os.path.exists(path):
        raise FileNotFoundError(f"no such file: {path}")
    name = str(path).lower()
    name = name.removesuffix(".gz")
    if name.endswith(".sdf"):
        rows = _read_sdf(path)
    elif name.endswith((".smi", ".ism")):
        rows = _read_smi(path)
    elif name.endswith((".csv", ".tsv", ".txt")):
        try:
            rows = _read_csv(path, smiles_column, id_column)
        except ValueError:
            if name.endswith(".txt"):
                rows = _read_smi(path)
            else:
                raise
    else:
        raise ValueError(f"unsupported file type: {Path(path).name} "
                         "(expected .csv, .tsv, .smi, .txt or .sdf, optionally .gz)")
    if limit:
        rows = rows[:limit]
    return rows


def load_mols(rows: list[dict], standardize_mols: bool = True) -> tuple[list[dict], list[dict]]:
    """Attach an RDKit molecule to each row. Returns (ok_rows, failed_rows).

    ok_rows carry ``mol`` (RDKit object, not serializable) and ``smiles`` as the canonical form.
    """
    ok, failed = [], []
    for row in rows:
        if row.get("error"):
            failed.append({"id": row["id"], "smiles": row.get("smiles", ""), "error": row["error"]})
            continue
        mol = parse(row["smiles"])
        if mol is None:
            failed.append({"id": row["id"], "smiles": row["smiles"], "error": "SMILES would not parse"})
            continue
        if standardize_mols:
            try:
                mol, _ = standardize(mol)
            except Exception as exc:  # noqa: BLE001 - best-effort; keep the original molecule
                row = dict(row, standardize_warning=str(exc))
        new = dict(row)
        new["mol"] = mol
        new["smiles"] = canonical(mol)
        ok.append(new)
    return ok, failed


def write_table(rows: list[dict], path: str, columns: list[str] | None = None) -> dict:
    """Write rows to CSV (or .sdf when the path says so). Returns a small summary for the caller."""
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        out.write_text("")
        return {"path": str(out), "rows": 0}

    if str(path).lower().endswith(".sdf"):
        writer = Chem.SDWriter(str(out))
        written = 0
        for row in rows:
            mol = row.get("mol") or parse(row.get("smiles", ""))
            if mol is None:
                continue
            mol = Chem.Mol(mol)
            mol.SetProp("_Name", str(row.get("id", "")))
            for key, value in row.items():
                if key in ("mol", "id") or value is None:
                    continue
                if isinstance(value, (dict, list)):
                    continue
                mol.SetProp(str(key), str(value))
            writer.write(mol)
            written += 1
        writer.close()
        return {"path": str(out), "rows": written}

    flat = []
    for row in rows:
        item = {k: v for k, v in row.items() if k != "mol" and not isinstance(v, (dict, list))}
        for key, value in (row.get("fields") or {}).items():
            item.setdefault(key, value)
        flat.append(item)
    cols = columns or list(dict.fromkeys(k for item in flat for k in item))
    with open(out, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(flat)
    return {"path": str(out), "rows": len(flat), "columns": cols}
