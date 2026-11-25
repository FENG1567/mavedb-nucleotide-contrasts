import re
import sys
import os
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from variant_map import fetch_cds, map_snv, parse_nt  # noqa: E402

ROOT = Path(os.environ.get('RL_PROJECT_ROOT',str(Path(__file__).resolve().parents[2]/'data/external/collection_project')))
WORK = ROOT / "work"
IN = ROOT / "incoming" / "multigene_scores"
OUT = WORK / "multigene_mapped"
OUT.mkdir(parents=True, exist_ok=True)

AA3 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q",
    "Glu": "E", "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K",
    "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W",
    "Tyr": "Y", "Val": "V", "Ter": "*",
}


def parse_pro(p):
    if not isinstance(p, str) or not p.startswith("p."):
        return None
    s = p.replace("p.(", "p.").replace(")", "")
    m = re.match(r"p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})$", s)
    if m:
        return int(m.group(2)), AA3.get(m.group(1)), AA3.get(m.group(3))
    m = re.match(r"p\.([A-Z])(\d+)([A-Z])$", s)
    if m:
        return int(m.group(2)), m.group(1), m.group(3)
    return None


def score_col(df):
    cols = list(df.columns)
    for c in ["score", "Score", "SCORE"]:
        if c in cols:
            return c
    cand = [c for c in cols if "score" in c.lower() and "error" not in c.lower()]
    return cand[0] if cand else None


def main():
    man = pd.read_csv(WORK / "multigene_manifest.csv", dtype={"label_prefix": str})
    summary = []
    for _, m in man.iterrows():
        urn = m.score_set
        fn = IN / (urn.replace(":", "_") + ".csv")
        if not fn.exists():
            summary.append({"score_set": urn, "gene": m.gene, "status": "missing"})
            continue
        head = pd.read_csv(fn, nrows=3, dtype=str, keep_default_na=False)
        sc = score_col(head)
        if sc is None:
            summary.append({"score_set": urn, "gene": m.gene, "status": "no_score_col"})
            continue
        df = pd.read_csv(fn, dtype=str, keep_default_na=False)
        df["_score"] = pd.to_numeric(df[sc], errors="coerce")
        df = df.dropna(subset=["_score"])
        rows = []
        kinds = Counter()
        for _, r in df.iterrows():
            pos = ref = alt = cls = None
            pro = r.get("hgvs_pro", "")
            if isinstance(pro, str) and pro and pro != "NA":
                pp = parse_pro(pro)
                if pp and pp[1] and pp[2]:
                    pos, ref, alt = pp
                    cls = "synonymous" if ref == alt else ("nonsense" if alt == "*" else "missense")
            if pos is None:
                nt = r.get("hgvs_nt", "")
                p = parse_nt(nt)
                if p is None:
                    kinds["missing"] += 1
                    continue
                kinds[p["kind"]] += 1
                if p["kind"] != "snv":
                    continue
                try:
                    cds = fetch_cds(p["acc"])
                    mm = map_snv(cds, p["pos"], p["ref"], p["alt"])
                except Exception:
                    kinds["cds_error"] += 1
                    continue
                if not mm or "error" in mm:
                    kinds["map_error"] += 1
                    continue
                pos, ref, alt, cls = mm["protein_pos"], mm["ref_aa"], mm["alt_aa"], mm["class"]
            rows.append({
                "score_set": urn, "gene": m.gene, "label_prefix": m.label_prefix,
                "variant_id": r.get("accession", ""),
                "hgvs_nt": r.get("hgvs_nt", ""), "hgvs_pro": r.get("hgvs_pro", ""),
                "score": r["_score"], "protein_pos": pos,
                "ref_aa": ref, "alt_aa": alt, "class": cls,
            })
        v = pd.DataFrame(rows)
        if len(v):
            v = v[v["class"] == "missense"].copy()
        safe = urn.replace(":", "_")
        v.to_csv(OUT / f"{safe}.csv", index=False)
        g = defaultdict(set)
        for _, r in v.iterrows():
            g[(int(r.protein_pos), r.ref_aa, r.alt_aa)].add(r.hgvs_nt or r.hgvs_pro or r.variant_id)
        n2 = sum(1 for x in g.values() if len(x) >= 2)
        summary.append({
            "score_set": urn, "gene": m.gene, "label_prefix": m.label_prefix,
            "status": "ok", "rows": len(df), "mapped": len(v),
            "groups": len(g), "groups_ge2": n2,
            "kinds": ";".join(f"{k}={v}" for k, v in kinds.items()),
        })
    s = pd.DataFrame(summary)
    s.to_csv(WORK / "multigene_map_summary.csv", index=False)
    print(s.groupby("status").size().to_string())
    print(s[s.status == "ok"].sort_values("groups_ge2", ascending=False).head(30).to_string(index=False))


if __name__ == "__main__":
    main()
