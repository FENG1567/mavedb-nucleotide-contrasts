"""Map MaveDB nucleotide HGVS variants to amino-acid changes.

Only exonic single-nucleotide substitutions are mapped in the main path.
Intronic, UTR, indel, and splice variants are labelled and excluded.
"""
from __future__ import annotations

import json
import os
import re
import urllib.request
from pathlib import Path

WORK = Path(os.environ.get("RL_WORK",str(Path(__file__).resolve().parents[2]/'data/external/collection_project/work')))
CDS_DIR = WORK / "cds"
CDS_DIR.mkdir(parents=True, exist_ok=True)

CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def translate(codon: str) -> str:
    return CODON_TABLE.get(codon.upper(), "?")


def fetch_cds(accession: str) -> str:
    """Fetch CDS for an Ensembl transcript or NCBI RefSeq mRNA, cached."""
    acc = accession.split(".")[0]
    cache = CDS_DIR / f"{acc}.txt"
    if cache.exists():
        return cache.read_text().strip()
    if acc.startswith("ENST"):
        url = f"https://rest.ensembl.org/sequence/id/{acc}?type=cds;content-type=text/plain"
        with urllib.request.urlopen(url, timeout=90) as r:
            seq = r.read().decode().strip()
    elif acc.startswith("NM_"):
        url = (
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
            f"?db=nuccore&id={acc}&rettype=fasta_cds_na&retmode=text"
        )
        with urllib.request.urlopen(url, timeout=90) as r:
            txt = r.read().decode()
        seq = "".join(line.strip() for line in txt.splitlines() if not line.startswith(">"))
    else:
        raise ValueError(f"unsupported accession {accession}")
    if not seq or set(seq.upper()) - set("ACGTN"):
        raise ValueError(f"bad CDS for {accession}: {seq[:50]}")
    cache.write_text(seq.upper() + "\n")
    return seq.upper()


NT_RE = re.compile(
    r"^(?P<acc>[^:]+):c\.(?P<pos>-?\d+|\*\d+)(?:(?P<off>[+-]\d+))?"
    r"(?:(?P<del>_(?P<pos2>-?\d+|\*\d+)))?"
    r"(?P<ref>[ACGTN])>(?P<alt>[ACGTN])$"
)
SNV_RE = re.compile(r"^c\.(?P<pos>\d+)(?P<ref>[ACGTN])>(?P<alt>[ACGTN])$")


def parse_nt(hgvs: str):
    """Return dict with kind and coordinates, or None if not a simple exonic SNV."""
    if not isinstance(hgvs, str) or not hgvs or hgvs == "NA":
        return None
    h = hgvs.strip()
    if ":c." not in h:
        return {"kind": "non_cds", "hgvs": h}
    acc, rest = h.split(":c.", 1)
    rest = "c." + rest
    # UTR
    if rest.startswith("c.-") or rest.startswith("c.*"):
        return {"kind": "utr", "hgvs": h}
    # intronic offset immediately after position
    m = re.match(r"^c\.(\d+)([+-])(\d+)([ACGTN])>([ACGTN])$", rest)
    if m:
        return {"kind": "splice", "hgvs": h, "pos": int(m.group(1)), "offset": m.group(2) + m.group(3)}
    m = SNV_RE.match(rest)
    if m:
        return {
            "kind": "snv",
            "hgvs": h,
            "acc": acc,
            "pos": int(m.group("pos")),
            "ref": m.group("ref"),
            "alt": m.group("alt"),
        }
    # any other exonic change (indel, delins, dup, etc.)
    if re.match(r"^c\.\d+(_\d+)?(del|ins|dup|delins)", rest):
        return {"kind": "indel", "hgvs": h}
    return {"kind": "other", "hgvs": h}


def map_snv(cds: str, pos: int, ref: str, alt: str):
    """Map a coding SNV (1-based CDS position) to codon/AA change."""
    if pos < 1 or pos > len(cds):
        return None
    if cds[pos - 1] != ref:
        return {"error": "ref_mismatch", "cds_ref": cds[pos - 1]}
    codon_start = ((pos - 1) // 3) * 3
    ref_codon = cds[codon_start:codon_start + 3]
    alt_codon = ref_codon[: (pos - 1) % 3] + alt + ref_codon[(pos - 1) % 3 + 1:]
    ref_aa = translate(ref_codon)
    alt_aa = translate(alt_codon)
    if ref_aa == alt_aa:
        cls = "synonymous"
    elif alt_aa == "*":
        cls = "nonsense"
    elif ref_aa == "*":
        cls = "stop_lost"
    else:
        cls = "missense"
    return {
        "cds_pos": pos,
        "codon_start": codon_start + 1,
        "ref_codon": ref_codon,
        "alt_codon": alt_codon,
        "ref_aa": ref_aa,
        "alt_aa": alt_aa,
        "protein_pos": codon_start // 3 + 1,
        "class": cls,
    }
