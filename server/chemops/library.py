"""Library-scale screening: the filter cascade, its attrition table, and deduplication.

The cascade reports what each stage removed and why. That attrition table is the part a chemist
actually argues with - a filter that throws away 80% of a deck is usually a mis-set threshold, not a
clean library, and you can only see that if the numbers are stage by stage.
"""

from __future__ import annotations

from rdkit import Chem

from . import alerts as alerts_mod
from . import mol as mol_mod
from . import props as props_mod

# Mutually exclusive outcomes recorded per compound.
PASS, FAIL_PROPERTY, FAIL_RULE, FAIL_ALERT = "pass", "property", "rule", "alert"


def dedupe(rows: list[dict]) -> tuple[list[dict], list[dict]]:
    """Collapse duplicates by InChIKey (the full key: stereochemistry counts as different).

    Returns (unique_rows, duplicate_records). The first occurrence wins and keeps a ``duplicate_ids``
    list so nothing is lost silently - knowing a compound appears five times in a deck matters.
    """
    seen: dict[str, dict] = {}
    duplicates: list[dict] = []
    for row in rows:
        key = mol_mod.inchikey(row["mol"]) or row["smiles"]
        if key in seen:
            seen[key].setdefault("duplicate_ids", []).append(row["id"])
            duplicates.append({"id": row["id"], "smiles": row["smiles"],
                               "duplicate_of": seen[key]["id"], "inchikey": key})
        else:
            seen[key] = row
    return list(seen.values()), duplicates


def _check_properties(d: dict, filters: dict) -> list[str]:
    """Return the list of property criteria this compound fails."""
    failures = []
    for key, bounds in filters.items():
        if key not in d or d[key] is None:
            continue  # descriptor unavailable (e.g. sa_score without the Contrib module)
        low, high = (list(bounds) + [None, None])[:2] if isinstance(bounds, (list, tuple)) else (None, bounds)
        value = d[key]
        if low is not None and value < low:
            failures.append(f"{key} {value} < {low}")
        if high is not None and value > high:
            failures.append(f"{key} {value} > {high}")
    return failures


def screen(rows: list[dict], filters: dict | None = None, require_rules: list[str] | None = None,
           alert_catalogs: list[str] | None = None, alert_action: str = "flag",
           deduplicate: bool = True, max_alerts: int = 0) -> dict:
    """Run the cascade over loaded rows.

    filters: ``{descriptor: [low, high]}``; either bound may be null. Keys are descriptor names from
      ``props.descriptors`` (mw, clogp, tpsa, hbd, hba, rotatable_bonds, qed, sa_score, ...).
    require_rules: rule-set names that must pass, e.g. ``["Lipinski Ro5", "Veber"]`` (prefix match).
    alert_catalogs: catalogs to evaluate (default PAINS, Brenk, NIH). ``alert_action`` is ``flag``
      (compute and report, keep the compound) or ``reject`` (fail compounds with more than
      ``max_alerts`` alerts).
    """
    filters = filters or {}
    require_rules = require_rules or []
    catalogs = alert_catalogs if alert_catalogs is not None else list(alerts_mod.DEFAULT_CATALOGS)
    if alert_action not in ("flag", "reject"):
        raise ValueError("alert_action must be 'flag' or 'reject'")

    stages: list[dict] = []
    n_in = len(rows)

    duplicates: list[dict] = []
    if deduplicate:
        rows, duplicates = dedupe(rows)
        stages.append({"stage": "deduplicate (InChIKey)", "in": n_in, "out": len(rows),
                       "removed": len(duplicates)})

    results: list[dict] = []
    reasons: dict[str, int] = {}
    for row in rows:
        m: Chem.Mol = row["mol"]
        d = props_mod.descriptors(m)
        rules = props_mod.rule_sets(m, d)
        alert_info = alerts_mod.alert_summary(m, catalogs) if catalogs else {"alert_count": 0,
                                                                             "alerts": [],
                                                                             "clean": True}
        prop_failures = _check_properties(d, filters)
        rule_failures = [
            r["rule"] for r in rules
            if not r["pass"] and any(r["rule"].lower().startswith(want.lower()) for want in require_rules)
        ]
        alert_reject = alert_action == "reject" and alert_info["alert_count"] > max_alerts

        if prop_failures:
            outcome, why = FAIL_PROPERTY, prop_failures
        elif rule_failures:
            outcome, why = FAIL_RULE, rule_failures
        elif alert_reject:
            outcome, why = FAIL_ALERT, [a["alert"] for a in alert_info["alerts"]]
        else:
            outcome, why = PASS, []

        for reason in why:
            reasons[reason] = reasons.get(reason, 0) + 1

        results.append({
            "id": row["id"],
            "smiles": row["smiles"],
            "outcome": outcome,
            "failed": why,
            "alert_count": alert_info["alert_count"],
            "alerts": [a["alert"] for a in alert_info["alerts"]],
            "rules_failed": [r["rule"] for r in rules if not r["pass"]],
            "mol": m,
            "fields": row.get("fields") or {},
            "duplicate_ids": row.get("duplicate_ids", []),
            **d,
        })

    passed = [r for r in results if r["outcome"] == PASS]
    after_props = [r for r in results if r["outcome"] not in (FAIL_PROPERTY,)]
    if filters:
        stages.append({"stage": f"property filters ({', '.join(filters)})", "in": len(rows),
                       "out": len(after_props), "removed": len(rows) - len(after_props)})
    after_rules = [r for r in after_props if r["outcome"] != FAIL_RULE]
    if require_rules:
        stages.append({"stage": f"rule sets ({', '.join(require_rules)})", "in": len(after_props),
                       "out": len(after_rules), "removed": len(after_props) - len(after_rules)})
    if alert_action == "reject":
        stages.append({"stage": f"structural alerts ({', '.join(catalogs)})", "in": len(after_rules),
                       "out": len(passed), "removed": len(after_rules) - len(passed)})

    flagged = sum(1 for r in results if r["alert_count"] > 0)
    return {
        "input_rows": n_in,
        "screened": len(rows),
        "passed": len(passed),
        "failed": len(results) - len(passed),
        "pass_rate": round(len(passed) / len(rows), 4) if rows else 0.0,
        "duplicates_removed": len(duplicates),
        "duplicates": duplicates[:50],
        "alert_flagged": flagged,
        "alert_action": alert_action,
        "catalogs": catalogs,
        "stages": stages,
        "top_reasons": sorted(reasons.items(), key=lambda kv: kv[1], reverse=True)[:15],
        "results": results,
    }


def summarize_distribution(results: list[dict], keys: tuple[str, ...] = (
        "mw", "clogp", "tpsa", "hbd", "hba", "rotatable_bonds", "qed", "sa_score")) -> dict:
    """Median and range for each descriptor across a set - the shape of the library in one glance."""
    out = {}
    for key in keys:
        values = sorted(r[key] for r in results if r.get(key) is not None)
        if not values:
            continue
        n = len(values)
        out[key] = {
            "min": values[0],
            "q1": values[n // 4],
            "median": values[n // 2],
            "q3": values[(3 * n) // 4],
            "max": values[-1],
            "n": n,
        }
    return out
