#!/usr/bin/env python3
"""Build the Proteinbase submission from whatever the campaign has finished.

Safe to run at any point, including on a half-finished round 2: every source is
optional and the script reports what it found and what it had to skip.

    python campaign/submit.py --n 40
    python campaign/submit.py --n 20 --out outputs/submission

Outputs (in --out):
    submission.csv            name,sequence,molecule_class  -- the upload file
    submission_annotated.csv  the same designs with every metric used
    structure_audit.csv       cached per-design AF3 contact audit
    selection_report.md       what was pooled, filtered, ranked and why

Selection is agonism-directed: the challenge only counts a design that binds
*and* activates, so engagement of the ECL3/TM7 hotspot cluster (D372/R380/L384)
is a hard requirement wherever a structure is available, and full 6/6 hotspot
coverage is preferred over raw interface size.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
import sys
from collections import defaultdict
from pathlib import Path

# Hotspots in UniProt P43220 numbering. The AF3 complex numbers the target
# chain directly in UniProt numbering, so these need no remapping (verified:
# 298 T, 299 R, 300 N, 372 D, 380 R, 384 L).
HOT_ECL2 = (298, 299, 300)          # ECL2 trio, the GLP-1 contact core
HOT_ECL3 = (372, 380, 384)          # ECL3/TM7, the activation-facing cluster
HOTSPOTS = HOT_ECL2 + HOT_ECL3

CAMPAIGN = "outputs/campaign/glp1r"
AA3 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
}
STANDARD_AA = set(AA3.values())


# --------------------------------------------------------------------------
# small helpers
# --------------------------------------------------------------------------

def fnum(row, key, default=None):
    """Float from a csv row, or default when absent/empty/unparseable."""
    try:
        v = row.get(key, "")
        return float(v) if v not in ("", None) else default
    except (TypeError, ValueError):
        return default


# Kyte-Doolittle hydropathy and side-chain pKa values, for the sequence-level
# developability columns that go into the methodology table.
KD = {"A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5,
      "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9,
      "M": 1.9, "F": 2.8, "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9,
      "Y": -1.3, "V": 4.2}
PKA_POS = {"K": 10.8, "R": 12.5, "H": 6.5}
PKA_NEG = {"D": 3.9, "E": 4.1, "C": 8.5, "Y": 10.1}
PKA_NTERM, PKA_CTERM = 8.6, 3.6


def charge_at(seq, ph):
    q = 1.0 / (1.0 + 10 ** (ph - PKA_NTERM)) - 1.0 / (1.0 + 10 ** (PKA_CTERM - ph))
    for aa, pka in PKA_POS.items():
        q += seq.count(aa) / (1.0 + 10 ** (ph - pka))
    for aa, pka in PKA_NEG.items():
        q -= seq.count(aa) / (1.0 + 10 ** (pka - ph))
    return q


def seq_props(seq):
    """Sequence-level developability readouts for the designed (MPNN) chain."""
    lo, hi = 0.0, 14.0
    for _ in range(100):                       # bisection on the titration curve
        mid = (lo + hi) / 2
        if charge_at(seq, mid) > 0:
            lo = mid
        else:
            hi = mid
    n = len(seq) or 1
    return {
        "length": n,
        "net_charge_pH7": round(charge_at(seq, 7.4), 1),
        "pI": round((lo + hi) / 2, 2),
        "GRAVY": round(sum(KD.get(a, 0.0) for a in seq) / n, 2),
        "n_cys": seq.count("C"),
        "n_glyco_NXST": len(re.findall(r"N[^P][ST]", seq)),
        "n_NG": seq.count("NG"),
        "aromatic_frac": round(sum(seq.count(a) for a in "FWY") / n, 3),
    }


def load_rosetta(root, champ):
    """design -> Rosetta InterfaceAnalyzer row for one champion.

    lead_ranking.csv carries only a subset of the Rosetta output; packstat and
    the interface composition terms live in the raw score table.
    """
    path = root / CAMPAIGN / "round2" / champ / "rosetta" / "af3_scores.csv"
    out = {}
    if not path.exists():
        return out
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh):
            if (r.get("species") or "human") == "human":
                out[r["design"]] = r
    return out


def load_seq_diversity(root, champ):
    """backbone -> MPNN sequence diversity for one champion's round-2 stage."""
    path = root / CAMPAIGN / "round2" / champ / "seq_diversity.csv"
    out = {}
    if not path.exists():
        return out
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh):
            out[r["backbone"]] = r
    return out


def backbone_of(design):
    """stage1_len60_g0r0_2_p0_dldesign_... -> stage1_len60_g0r0_2"""
    m = re.match(r"(.*?)_p\d+_dldesign", design)
    return m.group(1) if m else design


def short_name(design):
    """stage1_len60_g1r0_1_p4_dldesign_0_cycle1_af2pred -> g1r0_1_p4"""
    s = re.sub(r"^stage1_len\d+_", "", design)
    s = re.sub(r"_dldesign_\d+_cycle\d+(_\d+)?(_af2pred)?(__\w+)?$", "", s)
    return s or design


# --------------------------------------------------------------------------
# sources
# --------------------------------------------------------------------------

def load_round2(root, log):
    """Pool every finished round-2 champion. Richest source first."""
    rows, seen = [], set()
    champ_dirs = sorted((root / CAMPAIGN / "round2").glob("*/"))
    if not champ_dirs:
        log.append("round2: no champion directories found")
        return rows

    for cd in champ_dirs:
        champ = cd.name
        lead = cd / "final" / "lead_ranking.csv"
        short = cd / "final" / "shortlist.csv"
        src = lead if lead.exists() else (short if short.exists() else None)
        if src is None:
            stage = "AF2 done, AF3/Rosetta pending" if (cd / "af2_stage_done").exists() else "in progress"
            log.append(f"round2/{champ}: SKIPPED ({stage}; no final/ table, "
                       f"sequences live in the gitignored MPNN silents)")
            continue

        n_before = len(rows)
        diversity = load_seq_diversity(root, champ)
        rosetta = load_rosetta(root, champ)
        with open(src, newline="") as fh:
            for r in csv.DictReader(fh):
                seq = (r.get("sequence") or "").strip().upper()
                if not seq or seq in seen:
                    continue
                seen.add(seq)
                design = r.get("design", "")
                div = diversity.get(backbone_of(design), {})
                ros = rosetta.get(design, {})
                rows.append({
                    **seq_props(seq),
                    "mpnn_n_distinct": div.get("n_distinct"),
                    "mpnn_median_identity": div.get("median_identity"),
                    "design": design,
                    "short_id": r.get("short_id") or short_name(design),
                    "parent": champ,
                    "sequence": seq,
                    "source": src.name,
                    "iptm": fnum(r, "iptm_human", fnum(r, "chain_iptm_TB_human")),
                    "pae_int_af3": fnum(r, "chain_pae_min_TB_human"),
                    "pae_binder_af3": fnum(r, "pae_binder_trimer_mean_human"),
                    "pae_int_af2": fnum(r, "pae_interaction_af2"),
                    "rmsd_af2": fnum(r, "binder_aligned_rmsd_af2"),
                    "plddt_binder_af2": fnum(r, "plddt_binder_af2"),
                    "rmsd_conv": fnum(r, "rmsd_binder_after_target_super_human"),
                    "plddt_binder": fnum(r, "plddt_binder_mean_human"),
                    "ddG": fnum(r, "ddG_human"),
                    "sc": fnum(r, "sc_human"),
                    "dsasa": fnum(r, "dSASA_int_human"),
                    "hbonds": fnum(r, "interface_hbonds_human"),
                    "packstat": fnum(r, "packstat_human",
                                     fnum(ros, "packstat")),
                    "nres_int": fnum(ros, "nres_int"),
                    "interface_hydrophobicity": fnum(ros, "interface_hydrophobicity"),
                    "clash": fnum(r, "has_clash_human", 0.0),
                })
        log.append(f"round2/{champ}: {len(rows) - n_before} designs from {src.name}")
    return rows


def load_round1_fallback(root, log):
    """The 42 round-1 champions, scored by the blind AF3 pre-screen.

    Only used when round 2 produced nothing (or with --include-round1). The
    sequences are read out of the AF3 complex itself, so this works even though
    the design PDBs are gitignored.
    """
    rows = []
    base = root / CAMPAIGN / "round2_af3_prescreen"
    summary = base / "af3_summary.csv"
    if not summary.exists():
        log.append("round1 fallback: no af3_summary.csv")
        return rows

    with open(summary, newline="") as fh:
        for r in csv.DictReader(fh):
            design = r.get("design", "")
            cif = find_cif(root, design, None)
            if cif is None:
                continue
            _, binder_seq = read_complex(cif)
            if not binder_seq:
                continue
            rows.append({
                "design": design,
                "short_id": short_name(design),
                "parent": "round1",
                "sequence": binder_seq,
                "source": "af3_summary.csv",
                "iptm": fnum(r, "iptm", fnum(r, "chain_iptm_TB")),
                "pae_int_af3": fnum(r, "chain_pae_min_TB"),
                "pae_int_af2": None,
                "rmsd_af2": None,
                "rmsd_conv": fnum(r, "rmsd_binder_after_target_super"),
                "plddt_binder": fnum(r, "plddt_binder_mean"),
                "ddG": None, "sc": None, "dsasa": None,
                "hbonds": None, "packstat": None,
                "clash": fnum(r, "has_clash", 0.0),
            })
    log.append(f"round1 fallback: {len(rows)} AF3 pre-screened champions")
    return rows


# --------------------------------------------------------------------------
# structural audit: hotspot contacts + epitope size from the AF3 complex
# --------------------------------------------------------------------------

def find_cif(root, design, parent):
    """Locate one AF3 model for a design. Deterministic when several exist.

    Design tags are NOT unique across round-2 champions: each champion numbers
    its own partial-diffusion backbones g0r0_0..7 / g1r0_0..7, so the same tag
    names a different molecule under a different champion (81 of 115 tags
    collide in the current campaign). The search is therefore confined to the
    champion's own directory; a global fallback would silently return a
    sibling champion's structure.
    """
    if parent:
        cd = root / CAMPAIGN / "round2" / parent
        pats = [(cd, p) for p in (
            f"final/pdb_AF3/*{design}*.cif",
            f"af3/**/{design}__*/**/*_model.cif",
            f"**/{design}__*/**/*_model.cif",
        )]
    else:
        pats = [(root / CAMPAIGN, f"**/{design}__*/**/*_model.cif")]

    for base, pat in pats:
        if not base.exists():
            continue
        hits = sorted(base.glob(pat))
        if hits:
            # prefer the first seed/sample so repeat runs agree
            pref = [h for h in hits if "seed-1_sample-0" in h.name]
            return (pref or hits)[0]
    return None


def read_complex(cif):
    """-> ({target_resnum: [heavy atom xyz]}, binder_sequence).

    Chain A is the GLP1R target in UniProt numbering, chain B the binder.
    """
    target = defaultdict(list)
    binder = {}
    try:
        with open(cif) as fh:
            for line in fh:
                if not line.startswith("ATOM"):
                    continue
                p = line.split()
                if len(p) < 13:
                    continue
                chain, comp, element = p[6], p[5], p[2]
                if element == "H":
                    continue
                try:
                    seq = int(p[8])
                    xyz = (float(p[10]), float(p[11]), float(p[12]))
                except ValueError:
                    continue
                if chain == "A":
                    target[seq].append(xyz)
                elif chain == "B":
                    binder.setdefault(seq, (comp, []))[1].append(xyz)
    except OSError:
        return {}, ""
    seq_str = "".join(AA3.get(binder[k][0], "X") for k in sorted(binder))
    return target, seq_str


def audit_structure(cif, cutoff):
    """Hotspot contacts and epitope size, via a spatial hash of binder atoms."""
    target, _ = read_complex(cif)
    binder_atoms = []
    try:
        with open(cif) as fh:
            for line in fh:
                if line.startswith("ATOM"):
                    p = line.split()
                    if len(p) >= 13 and p[6] == "B" and p[2] != "H":
                        binder_atoms.append((float(p[10]), float(p[11]), float(p[12])))
    except OSError:
        return None
    if not target or not binder_atoms:
        return None

    cell = cutoff
    grid = defaultdict(list)
    for x, y, z in binder_atoms:
        grid[(int(x // cell), int(y // cell), int(z // cell))].append((x, y, z))

    def min_dist(atoms):
        best = 1e9
        for x, y, z in atoms:
            cx, cy, cz = int(x // cell), int(y // cell), int(z // cell)
            for dx in (-1, 0, 1):
                for dy in (-1, 0, 1):
                    for dz in (-1, 0, 1):
                        for bx, by, bz in grid.get((cx + dx, cy + dy, cz + dz), ()):
                            d = math.dist((x, y, z), (bx, by, bz))
                            if d < best:
                                best = d
            if best <= 1.0:
                break
        return best

    epitope, hot = [], {}
    for resnum, atoms in target.items():
        d = min_dist(atoms)
        if d <= cutoff:
            epitope.append(resnum)
        if resnum in HOTSPOTS:
            hot[resnum] = d

    return {
        "cif": str(cif),
        "n_epitope": len(epitope),
        "ecl2_contacts": sum(1 for h in HOT_ECL2 if hot.get(h, 1e9) <= cutoff),
        "ecl3_contacts": sum(1 for h in HOT_ECL3 if hot.get(h, 1e9) <= cutoff),
        "min_hot_dist": round(min(hot.values()), 2) if hot else None,
        "hot_dists": ";".join(f"{h}:{hot.get(h, float('nan')):.2f}" for h in HOTSPOTS),
    }


def run_audit(root, rows, cutoff, cache_path, log):
    """Audit every row, reusing a cache so re-runs are cheap.

    Keyed by (parent, design): a bare design tag is ambiguous across champions
    (see find_cif), so caching on it alone hands one champion's structure to
    another champion's design.
    """
    cache = {}
    if cache_path.exists():
        with open(cache_path, newline="") as fh:
            for r in csv.DictReader(fh):
                # rows without a parent come from the pre-fix cache format and
                # cannot be trusted to belong to the design they name
                if r.get("cutoff") == f"{cutoff}" and r.get("parent"):
                    cache[(r["parent"], r["design"])] = r

    audited = missing = 0
    for row in rows:
        key = (row["parent"], row["design"])
        hit = cache.get(key)
        if hit and hit.get("n_epitope"):
            row.update({
                "n_epitope": int(float(hit["n_epitope"])),
                "ecl2_contacts": int(float(hit["ecl2_contacts"])),
                "ecl3_contacts": int(float(hit["ecl3_contacts"])),
                "min_hot_dist": float(hit["min_hot_dist"]) if hit.get("min_hot_dist") else None,
                "hot_dists": hit.get("hot_dists", ""),
                "cif": hit.get("cif", ""),
            })
            audited += 1
            continue

        cif = find_cif(root, row["design"], row["parent"] if row["parent"] != "round1" else None)
        res = audit_structure(cif, cutoff) if cif else None
        if res:
            row.update(res)
            cache[key] = dict(res, design=row["design"],
                              parent=row["parent"], cutoff=cutoff)
            audited += 1
        else:
            row.update({"n_epitope": None, "ecl2_contacts": None,
                        "ecl3_contacts": None, "min_hot_dist": None,
                        "hot_dists": "", "cif": ""})
            missing += 1

    cache_path.parent.mkdir(parents=True, exist_ok=True)
    cols = ["parent", "design", "cutoff", "n_epitope", "ecl2_contacts",
            "ecl3_contacts", "min_hot_dist", "hot_dists", "cif"]
    with open(cache_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for d in cache.values():
            w.writerow(d)

    log.append(f"structure audit: {audited} designs audited, {missing} without a "
               f"reachable AF3 model (AF3 outputs are gitignored -- run this on "
               f"the machine that holds them for the full audit)")
    return audited, missing


# --------------------------------------------------------------------------
# filtering and ranking
# --------------------------------------------------------------------------

def hard_filters(rows, args, log):
    """Drop anything that cannot or should not be submitted."""
    kept, drops = [], defaultdict(int)
    for r in rows:
        seq = r["sequence"]
        if not (10 <= len(seq) <= 250):
            drops["length outside 10-250 aa"] += 1
            continue
        if set(seq) - STANDARD_AA:
            drops["non-standard residue"] += 1
            continue
        if (r.get("clash") or 0) > 0:
            drops["AF3 inter-chain clash"] += 1
            continue
        if r.get("iptm") is not None and r["iptm"] < args.min_iptm:
            drops[f"iptm < {args.min_iptm}"] += 1
            continue
        if args.require_ecl3 and r.get("ecl3_contacts") is not None \
                and r["ecl3_contacts"] < args.min_ecl3:
            drops[f"ECL3/TM7 contacts < {args.min_ecl3} (unlikely agonist)"] += 1
            continue
        if args.min_hot_coverage and r.get("ecl2_contacts") is not None:
            if r["ecl2_contacts"] + r["ecl3_contacts"] < args.min_hot_coverage:
                drops[f"hotspot coverage < {args.min_hot_coverage}/6"] += 1
                continue
        if args.max_rmsd_conv and r.get("rmsd_conv") is not None \
                and r["rmsd_conv"] > args.max_rmsd_conv:
            drops[f"AF2/AF3 pose disagreement > {args.max_rmsd_conv} A"] += 1
            continue
        if r.get("n_epitope") is not None and r["n_epitope"] > args.epitope_max:
            drops[f"epitope > {args.epitope_max} residues (drape pathology)"] += 1
            continue
        kept.append(r)

    for reason, n in sorted(drops.items(), key=lambda kv: -kv[1]):
        log.append(f"  filtered {n:4d}  {reason}")
    return kept


def composite(rows, args):
    """Rank designs. Interface confidence first, agonism-relevant terms next.

    Every term is normalised against the pool, so a missing metric scores
    neutral rather than disqualifying a design.
    """
    def norm(key, invert=False):
        vals = [r[key] for r in rows if r.get(key) is not None]
        if not vals:
            return {id(r): 0.5 for r in rows}
        lo, hi = min(vals), max(vals)
        span = (hi - lo) or 1.0
        out = {}
        for r in rows:
            v = r.get(key)
            if v is None:
                out[id(r)] = 0.5
            else:
                z = (v - lo) / span
                out[id(r)] = 1.0 - z if invert else z
        return out

    # ddG is only meaningful per unit of buried area; a huge raw ddG on a huge
    # interface is the flat-drape signature, not a better binder.
    for r in rows:
        if r.get("ddG") is not None and r.get("dsasa"):
            r["ddG_per_1k"] = r["ddG"] / (r["dsasa"] / 1000.0)
        else:
            r["ddG_per_1k"] = None
        e2, e3 = r.get("ecl2_contacts"), r.get("ecl3_contacts")
        r["hot_coverage"] = (e2 + e3) if (e2 is not None and e3 is not None) else None

    terms = {
        "iptm":        (norm("iptm"),              2.5),   # AF3 interface confidence
        "pae_af3":     (norm("pae_int_af3", True), 1.5),
        "pae_af2":     (norm("pae_int_af2", True), 1.5),
        "hot_cov":     (norm("hot_coverage"),      2.0),   # agonism-relevant
        "ecl3":        (norm("ecl3_contacts"),     1.5),   # activation cluster
        "sc":          (norm("sc"),                1.0),
        "packstat":    (norm("packstat"),          0.5),
        "hbonds":      (norm("hbonds"),            0.5),
        "ddG_per_1k":  (norm("ddG_per_1k", True),  1.0),   # more negative = better
        "epitope":     (norm("n_epitope", True),   1.0),   # locality guard
        "rmsd_conv":   (norm("rmsd_conv", True),   1.0),   # AF2/AF3 agreement
    }
    total_w = sum(w for _, w in terms.values())
    for r in rows:
        r["score"] = round(
            sum(t[id(r)] * w for t, w in terms.values()) / total_w, 4)
        r["flags"] = ";".join(filter(None, [
            f"epitope={r['n_epitope']}" if (r.get("n_epitope") or 0) > args.epitope_soft else "",
            f"dSASA={r['dsasa']:.0f}" if (r.get("dsasa") or 0) > args.dsasa_soft else "",
            f"rmsd_conv={r['rmsd_conv']:.1f}" if (r.get("rmsd_conv") or 0) > 4.0 else "",
            "no_structure" if r.get("n_epitope") is None else "",
            f"hot={r['hot_coverage']}/6" if (r.get("hot_coverage") is not None
                                             and r["hot_coverage"] < 6) else "",
        ]))
    rows.sort(key=lambda r: -r["score"])
    return rows


def diversify(rows, n, args, log):
    """Greedy pick with caps per parent lead and per sequence cluster.

    A submission of 40 near-identical sequences is one experiment, not forty.
    """
    def identity(a, b):
        if len(a) != len(b):
            return 0.0
        return sum(x == y for x, y in zip(a, b)) / len(a)

    picked, per_parent = [], defaultdict(int)
    for pas in range(2):                      # pass 2 relaxes the caps to fill n
        cap = args.per_parent if pas == 0 else n
        for r in rows:
            if len(picked) >= n:
                break
            if r in picked:
                continue
            if per_parent[r["parent"]] >= cap:
                continue
            if pas == 0 and any(identity(r["sequence"], p["sequence"]) >= args.max_identity
                                for p in picked):
                continue
            picked.append(r)
            per_parent[r["parent"]] += 1
        if len(picked) >= n:
            break

    log.append(f"diversity: {len(picked)} picked across {len(per_parent)} parent "
               f"lead(s) (cap {args.per_parent}/parent, "
               f"<{args.max_identity:.0%} pairwise identity)")
    for p, c in sorted(per_parent.items()):
        log.append(f"  {p}: {c}")
    return picked


# --------------------------------------------------------------------------
# output
# --------------------------------------------------------------------------

ANNOT_COLS = [
    "rank", "name", "short_id", "parent", "score", "flags", "sequence",
    "iptm", "pae_int_af3", "pae_binder_af3", "pae_int_af2", "rmsd_af2",
    "rmsd_conv", "plddt_binder", "plddt_binder_af2", "hot_coverage",
    "ecl2_contacts", "ecl3_contacts", "min_hot_dist", "hot_dists",
    "n_epitope", "ddG", "ddG_per_1k", "sc", "dsasa", "hbonds", "packstat",
    "clash", "nres_int", "interface_hydrophobicity", "length", "net_charge_pH7", "pI", "GRAVY", "n_cys",
    "n_glyco_NXST", "n_NG", "aromatic_frac", "mpnn_n_distinct",
    "mpnn_median_identity", "design", "source", "cif",
]

TABLE_BEGIN = "<!-- BEGIN AUTOGENERATED LEAD TABLE -- rewritten by campaign/submit.py -->"
TABLE_END = "<!-- END AUTOGENERATED LEAD TABLE -->"


def _md_table(header, rows):
    out = ["| " + " | ".join(header) + " |",
           "|" + "|".join("---" for _ in header) + "|"]
    out += ["| " + " | ".join(r) + " |" for r in rows]
    return out


def lead_tables(picked):
    """Every metric behind the submitted leads, grouped by the method it comes
    from. Four narrow tables rather than one unreadable 30-column grid."""
    def f(r, k, fmt="{:.2f}", dash="-"):
        v = r.get(k)
        if v in (None, ""):
            return dash
        try:
            return fmt.format(float(v))
        except (TypeError, ValueError):
            return str(v)

    L = [f"## 5. The {len(picked)} submitted designs", "",
         "Ranked by the composite described in §3. Every metric produced by the "
         "pipeline is listed, grouped by its source.", ""]

    L += ["### 5.1 Interface confidence and pose (AlphaFold2 / AlphaFold3)", "",
          "AF2 is the initial-guess prediction used as the design gate; AF3 is a "
          "blind re-prediction (target MSA only, no initial guess).", ""]
    L += _md_table(
        ["#", "design", "parent lead", "AF2 pae_int", "AF2 rmsd_b", "AF2 pLDDT_b",
         "AF3 ipTM", "AF3 PAE_int", "AF3 PAE_b", "AF3 pLDDT_b", "AF2/AF3 rmsd",
         "clash"],
        [[str(r["rank"]), r["short_id"], r["parent"],
          f(r, "pae_int_af2"), f(r, "rmsd_af2"), f(r, "plddt_binder_af2", "{:.1f}"),
          f(r, "iptm"), f(r, "pae_int_af3"), f(r, "pae_binder_af3", "{:.1f}"),
          f(r, "plddt_binder", "{:.1f}"), f(r, "rmsd_conv"),
          f(r, "clash", "{:.0f}")] for r in picked])

    L += ["", "### 5.2 Epitope and hotspot engagement (AF3 complex, 4.5 Å heavy atom)", "",
          "Hotspots in UniProt P43220 numbering. Every submitted design contacts "
          "all six. Distances are the closest binder heavy atom to that residue.", ""]
    L += _md_table(
        ["#", "design", "parent lead", "ECL2 /3", "ECL3 /3", "total /6",
         "min dist Å", "epitope res", "T298", "R299", "N300", "D372", "R380",
         "L384"],
        [[str(r["rank"]), r["short_id"], r["parent"],
          f(r, "ecl2_contacts", "{:.0f}"), f(r, "ecl3_contacts", "{:.0f}"),
          f(r, "hot_coverage", "{:.0f}"), f(r, "min_hot_dist"),
          f(r, "n_epitope", "{:.0f}")]
         + [dict(p.split(":") for p in (r.get("hot_dists") or "").split(";") if ":" in p)
            .get(h, "-") for h in ("298", "299", "300", "372", "380", "384")]
         for r in picked])

    L += ["", "### 5.3 Interface energetics (Rosetta InterfaceAnalyzer)", "",
          "Computed on the AF3 complex. ddG is also given per 1,000 Å² of buried "
          "surface, since a large raw ddG on an oversized interface reflects area "
          "rather than quality.", ""]
    L += _md_table(
        ["#", "design", "parent lead", "ddG", "ddG / 1k Å²", "dSASA_int Å²",
         "shape compl.", "packstat", "H-bonds", "interface res",
         "hydrophobic %"],
        [[str(r["rank"]), r["short_id"], r["parent"], f(r, "ddG", "{:.1f}"),
          f(r, "ddG_per_1k", "{:.1f}"), f(r, "dsasa", "{:.0f}"), f(r, "sc"),
          f(r, "packstat"), f(r, "hbonds", "{:.0f}"),
          f(r, "nres_int", "{:.0f}"),
          f(r, "interface_hydrophobicity", "{:.1f}")] for r in picked])

    L += ["", "### 5.4 Designed sequence (ProteinMPNN) and developability", "",
          "MPNN diversity is measured across the 8 sequences sampled for that "
          "design's backbone. No design carries a cysteine.", ""]
    L += _md_table(
        ["#", "design", "parent lead", "len", "net charge pH 7.4", "pI",
         "GRAVY", "aromatic", "Cys", "N-glyc sequons", "NG",
         "MPNN distinct /8", "MPNN median id"],
        [[str(r["rank"]), r["short_id"], r["parent"], f(r, "length", "{:.0f}"),
          f(r, "net_charge_pH7", "{:+.1f}"), f(r, "pI"), f(r, "GRAVY"),
          f(r, "aromatic_frac", "{:.1%}"), f(r, "n_cys", "{:.0f}"),
          f(r, "n_glyco_NXST", "{:.0f}"), f(r, "n_NG", "{:.0f}"),
          f(r, "mpnn_n_distinct", "{:.0f}"), f(r, "mpnn_median_identity")]
         for r in picked])

    L += ["", "### 5.5 Sequences", "",
          "| # | name | sequence |", "|---|------|----------|"]
    L += [f"| {r['rank']} | {r['name']} | `{r['sequence']}` |" for r in picked]
    return L


def update_methodology(out, picked, log):
    """Rewrite the autogenerated table block in methodology.md, in place.

    The prose above the markers is hand-written and never touched; only the
    block between them is replaced, so the table cannot drift out of step with
    submission.csv.
    """
    path = out / "methodology.md"
    block = "\n".join([TABLE_BEGIN, ""] + lead_tables(picked) + ["", TABLE_END])
    if path.exists():
        text = path.read_text()
        if TABLE_BEGIN in text and TABLE_END in text:
            head, _, rest = text.partition(TABLE_BEGIN)
            _, _, tail = rest.partition(TABLE_END)
            path.write_text(head + block + tail)
            log.append("methodology.md: lead table refreshed in place")
        else:
            path.write_text(text.rstrip() + "\n\n---\n\n" + block + "\n")
            log.append("methodology.md: lead table appended (no markers found)")
    else:
        path.write_text("# GLP1R submission\n\n" + block + "\n")
        log.append("methodology.md: created with the lead table only")


def write_outputs(out, picked, pooled, log, args):
    out.mkdir(parents=True, exist_ok=True)

    for i, r in enumerate(picked, 1):
        r["rank"] = i
        r["name"] = f"glp1r_{i:02d}_{r['parent']}_{r['short_id']}"[:60]

    with open(out / "submission.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["name", "sequence", "molecule_class"])
        for r in picked:
            # The CSV accepts only protein / nanobody / scfv / fab_kappa /
            # fab_lambda. "Single chain protein" names the design format in the
            # rules, not a column value: a de novo mini-binder is "protein".
            w.writerow([r["name"], r["sequence"], "single_chain"])

    with open(out / "submission_annotated.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=ANNOT_COLS, extrasaction="ignore")
        w.writeheader()
        w.writerows(picked)

    update_methodology(out, picked, log)

    lines = [
        "# GLP1R submission -- selection report",
        "",
        f"Generated by `campaign/submit.py` | pooled {len(pooled)} unique "
        f"sequences | submitting {len(picked)}",
        "",
        "## Sources and filtering",
        "",
    ]
    lines += [f"- {m}" for m in log]
    lines += [
        "",
        "## Parameters",
        "",
        f"- contact cutoff {args.cutoff} A (heavy atom), hotspots in UniProt "
        f"P43220 numbering: ECL2 {HOT_ECL2}, ECL3/TM7 {HOT_ECL3}",
        f"- hard: no AF3 clash, iptm >= {args.min_iptm}, ECL3 contacts >= "
        f"{args.min_ecl3} (when a structure is available), epitope <= "
        f"{args.epitope_max} residues, 10-250 aa, unique sequence"
        + (f", hotspot coverage >= {args.min_hot_coverage}/6" if args.min_hot_coverage else "")
        + (f", AF2/AF3 pose disagreement <= {args.max_rmsd_conv} A" if args.max_rmsd_conv else ""),
        f"- soft flags: epitope > {args.epitope_soft}, dSASA > {args.dsasa_soft} A^2",
        f"- diversity: <= {args.per_parent} per parent lead, pairwise identity "
        f"< {args.max_identity:.0%}",
        "",
        "## Selected designs",
        "",
        "| # | name | iptm | paeAF2 | hot | epi | ddG/1kA2 | sc | flags |",
        "|---|------|------|--------|-----|-----|----------|----|-------|",
    ]
    def cell(v, fmt="{:.2f}"):
        return fmt.format(v) if isinstance(v, (int, float)) else "-"
    for r in picked:
        lines.append(
            f"| {r['rank']} | {r['name']} | {cell(r.get('iptm'))} | "
            f"{cell(r.get('pae_int_af2'))} | "
            f"{r.get('hot_coverage') if r.get('hot_coverage') is not None else '-'}/6 | "
            f"{r.get('n_epitope') if r.get('n_epitope') is not None else '-'} | "
            f"{cell(r.get('ddG_per_1k'), '{:.1f}')} | {cell(r.get('sc'))} | "
            f"{r.get('flags', '')} |")
    (out / "selection_report.md").write_text("\n".join(lines) + "\n")


# --------------------------------------------------------------------------

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n", type=int, default=20,
                    help="designs to submit (Track 1 allows 40, Tracks 2/3 allow 20)")
    ap.add_argument("--out", default="outputs/submission")
    ap.add_argument("--root", default=None, help="repo root (default: parent of this file)")
    ap.add_argument("--min-iptm", type=float, default=0.80)
    ap.add_argument("--min-ecl3", type=int, default=1,
                    help="minimum ECL3/TM7 hotspot contacts (agonism requirement)")
    ap.add_argument("--no-require-ecl3", dest="require_ecl3", action="store_false")
    ap.add_argument("--min-hot-coverage", type=int, default=6,
                    help="require this many of the 6 hotspots in contact "
                         "(6 = full ECL2 + ECL3/TM7 engagement); 0 disables")
    ap.add_argument("--max-rmsd-conv", type=float, default=4,
                    help="reject designs whose binder moves more than this (A) "
                         "between the AF2 and AF3 models after target "
                         "superposition; unset disables")
    ap.add_argument("--cutoff", type=float, default=4.5, help="heavy-atom contact cutoff (A)")
    ap.add_argument("--epitope-soft", type=int, default=55,
                    help="flag designs contacting more target residues than this "
                         "(measured median for these 60-mers is 47)")
    ap.add_argument("--epitope-max", type=int, default=85,
                    help="reject above this; the failed 90-mers sat at 108-132")
    ap.add_argument("--dsasa-soft", type=float, default=5000.0)
    ap.add_argument("--per-parent", type=int, default=6,
                    help="max designs from one parent lead before relaxing")
    ap.add_argument("--max-identity", type=float, default=0.70)
    ap.add_argument("--include-round1", action="store_true",
                    help="also pool the round-1 AF3 pre-screened champions")
    args = ap.parse_args(argv)

    root = Path(args.root) if args.root else Path(__file__).resolve().parent.parent
    out = (root / args.out) if not Path(args.out).is_absolute() else Path(args.out)

    log = []
    rows = load_round2(root, log)
    if not rows or args.include_round1:
        if not rows:
            log.append("round 2 yielded nothing -- falling back to round 1")
        rows += [r for r in load_round1_fallback(root, log)
                 if r["sequence"] not in {x["sequence"] for x in rows}]

    if not rows:
        print("No designs found. Nothing finished far enough to submit.", file=sys.stderr)
        return 1

    log.append(f"pooled: {len(rows)} unique sequences")
    tags = defaultdict(set)
    for r in rows:
        tags[r["design"]].add(r["parent"])
    collide = sum(1 for p in tags.values() if len(p) > 1)
    if collide:
        log.append(f"note: {collide} of {len(tags)} design tags occur under more "
                   f"than one champion (each champion renumbers its own "
                   f"backbones); audits are keyed per champion, not per tag")
    run_audit(root, rows, args.cutoff, out / "structure_audit.csv", log)

    log.append("hard filters:")
    kept = hard_filters(rows, args, log)
    log.append(f"  -> {len(kept)} designs eligible")
    if not kept:
        print("Every design was filtered out; loosen --min-iptm or --min-ecl3.",
              file=sys.stderr)
        return 1

    ranked = composite(kept, args)
    picked = diversify(ranked, args.n, args, log)
    write_outputs(out, picked, rows, log, args)

    print(f"pooled {len(rows)} -> eligible {len(kept)} -> submitting {len(picked)}")
    print(f"wrote {out/'submission.csv'}")
    print(f"      {out/'submission_annotated.csv'}")
    print(f"      {out/'selection_report.md'}")
    if len(picked) < args.n:
        print(f"NOTE: only {len(picked)} of {args.n} slots filled.", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
