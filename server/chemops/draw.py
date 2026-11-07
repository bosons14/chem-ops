"""2D depiction to PNG or SVG, because a chemist judges a hit list by looking at the structures."""

from __future__ import annotations

from pathlib import Path

from rdkit import Chem
from rdkit.Chem import Draw, rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D


def _prepared(mol: Chem.Mol) -> Chem.Mol:
    mol = Chem.Mol(mol)
    rdDepictor.Compute2DCoords(mol)
    rdDepictor.StraightenDepiction(mol)
    return mol


def grid(rows: list[dict], path: str, per_row: int = 4, size: int = 300,
         legend_keys: tuple[str, ...] = ("id",), highlight_key: str | None = None,
         align_to: str | None = None) -> dict:
    """Draw a grid of structures with legends.

    legend_keys: row fields to put under each structure, e.g. ``("id", "mw", "qed")``.
    highlight_key: row field holding atom indices to highlight (alert atoms, for instance).
    align_to: SMARTS or SMILES core; when given, every molecule is oriented on that core so a series
      can be compared visually, which is the whole point of looking at a series.
    """
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise ValueError("nothing to draw")

    mols, legends, highlights = [], [], []
    template = None
    if align_to:
        template = Chem.MolFromSmarts(align_to) or Chem.MolFromSmiles(align_to)
        if template is not None:
            rdDepictor.Compute2DCoords(template)

    for row in rows:
        m = row.get("mol") or Chem.MolFromSmiles(row.get("smiles", ""))
        if m is None:
            continue
        m = _prepared(m)
        if template is not None and m.HasSubstructMatch(template):
            rdDepictor.GenerateDepictionMatching2DStructure(m, template, acceptFailure=True)
        mols.append(m)
        parts = []
        for key in legend_keys:
            value = row.get(key)
            if value not in (None, "", []):
                parts.append(str(value) if key == "id" else f"{key}={value}")
        legends.append(" | ".join(parts))
        highlights.append([int(i) for i in (row.get(highlight_key) or [])] if highlight_key else [])

    svg = str(out).lower().endswith(".svg")
    img = Draw.MolsToGridImage(
        mols, molsPerRow=min(per_row, len(mols)), subImgSize=(size, size),
        legends=legends, highlightAtomLists=highlights if highlight_key else None,
        useSVG=svg,
    )
    data = img.data if hasattr(img, "data") else img
    if svg:
        out.write_text(data if isinstance(data, str) else data.decode())
    else:
        payload = data if isinstance(data, bytes) else getattr(img, "data", None)
        if payload is None:  # PIL image
            img.save(str(out))
        else:
            out.write_bytes(payload)
    return {"path": str(out), "molecules": len(mols), "aligned": bool(template is not None)}


def single(mol: Chem.Mol, path: str, size: int = 400, highlight_atoms: list[int] | None = None,
           legend: str = "") -> dict:
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    m = _prepared(mol)
    svg = str(out).lower().endswith(".svg")
    drawer = (rdMolDraw2D.MolDraw2DSVG if svg else rdMolDraw2D.MolDraw2DCairo)(size, size)
    rdMolDraw2D.PrepareAndDrawMolecule(
        drawer, m, legend=legend,
        highlightAtoms=[int(i) for i in (highlight_atoms or [])] or None,
    )
    drawer.FinishDrawing()
    text = drawer.GetDrawingText()
    out.write_text(text) if svg else out.write_bytes(text)
    return {"path": str(out)}
