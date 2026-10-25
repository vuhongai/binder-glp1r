#!/usr/bin/env python3
"""AF3 pre-screen of the round-1 champions, BEFORE committing to round 2.

The round-1 funnel only applied the AF2 initial-guess gate (pae/rmsd/plddt).
AF3 is the blind, independent check: it folds the champion's binder sequence
against the target with NO initial guess, so a champion that holds up here
(high iptm, low chain pae, high binder plddt) is a much safer basin to spend
round-2 diffusion on.

This screens the champion PDBs (chain A = binder sequence, target from the
FASTA) and reports a per-champion pass/fail against the round-2 AF3
thresholds.  It uses the SAME target MSA cache round 2 will use
(results/<target>/af3_msa_cache), so running it first costs round 2 nothing
extra.

    python campaign/af3_prescreen.py campaign/glp1r_round2/round2.yaml \
        --champions-file campaign/champions_round2.txt --gpu 1 --check
    python campaign/af3_prescreen.py campaign/glp1r_round2/round2.yaml \
        --champions-file campaign/champions_round2.txt --gpu 1

Needs a FREE GPU (AF3 preallocates 95% of a card); the LLM server must be
stopped on that card first.
"""
from __future__ import annotations

import argparse
import csv
import os

import _campaign as C


def _champions(path: str) -> list[tuple[str, str]]:
    """tag=source_pdb pairs, one per line (the champions_round2.txt format)."""
    out = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            tag, _, src = line.partition("=")
            src = src if os.path.isabs(src) else os.path.join(C.ROOT, src)
            out.append((tag, src))
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("config")
    ap.add_argument("--champions-file", required=True,
                    help="tag=source_pdb pairs, one per line")
    ap.add_argument("--gpu", type=int, default=None,
                    help="GPU to run on (default: first configured)")
    ap.add_argument("--out", default=None,
                    help="output dir (default <run_root>/round2_af3_prescreen)")
    ap.add_argument("--check", action="store_true")
    args = ap.parse_args(argv)

    cfg = C.load(args.config)
    C.require_inputs(cfg)
    from binder.tools import af3 as af3tool
    champions = _champions(args.champions_file)
    if not champions:
        raise SystemExit(f"no champions in {args.champions_file}")
    pdbs = [src for _tag, src in champions]
    missing = [p for p in pdbs if not os.path.isfile(p)]
    if missing:
        raise SystemExit(f"champion PDBs missing:\n  " + "\n  ".join(missing))

    gpu = args.gpu if args.gpu is not None else cfg.resources.gpus[0]
    out_root = args.out or os.path.join(cfg.run_root, "round2_af3_prescreen")
    msa_cache = os.path.join(cfg.run_root, "af3_msa_cache")
    species = list(cfg.target.species)

    print(f"champions     {len(champions)}")
    print(f"species       {species}")
    print(f"gpu           {gpu}")
    print(f"msa_cache     {msa_cache} (shared with round 2)")
    print(f"out           {out_root}")
    if args.check:
        print("check only -- not running AF3")
        return 0

    log = os.path.join(out_root, "prescreen.log")
    notifier = C.setup(cfg, log=log)
    with notifier.arm(f"{cfg.target.name} af3_prescreen", log):
        af3 = af3tool.AF3(cfg)
        af3.check()
        # 1. target MSA, per species, resumable (shared with round 2).
        for i, sp in enumerate(species):
            cache_file = os.path.join(msa_cache, f"{sp}_target_data.json")
            if os.path.exists(cache_file):
                print(f"AF3 MSA cache present for {sp}: {cache_file}")
            else:
                print(f"AF3 MSA precompute for {sp} (GPU {gpu})")
                af3.precompute_target_msa(sp, gpu, msa_cache)
        # 2. inputs: one JSON per (champion, species).
        af3.prepare_inputs(pdbs, out_root, species_list=species,
                           msa_cache_dir=msa_cache, binder_chain="A")
        # 3. screen: one species per GPU (here a single species, one GPU).
        for sp in species:
            af3.screen(sp, gpu, out_root)
        # 4. summary against the round-2 AF3 thresholds (no rmsd join: the
        #    blind pose is independent of the champion pose).
        thr = cfg.thresholds.af3
        rows = af3.summary(out_root,
                           out_csv=os.path.join(out_root, "af3_summary.csv"),
                           species_list=species, rmsd_max=thr.rmsd_max,
                           pass_list_out=os.path.join(out_root,
                                                      "af3_pass.txt"))

    # Per-champion table, best-first, for the decision.
    rows.sort(key=lambda r: (float(r.get("iptm") or 0),
                             float(r.get("chain_pae_min_TB") or 999)))
    print(f"\n{'champion':16s} {'iptm':>6s} {'pae_TB':>7s} "
          f"{'plddt_b':>8s}  pass")
    passed = []
    for r in rows:
        ok = os.path.exists(os.path.join(out_root, "af3_pass.txt")) and any(
            line.split(",")[0] == r["design"]
            for line in open(os.path.join(out_root, "af3_pass.txt"))
            if line.strip())
        flag = "PASS" if ok else "fail"
        if ok:
            passed.append(r["design"])
        print(f"{r['design']:16s} {r.get('iptm','?'):>6s} "
              f"{r.get('chain_pae_min_TB','?'):>7s} "
              f"{r.get('plddt_binder_mean','?'):>8s}  {flag}")
    print(f"\n{len(passed)}/{len(rows)} champions pass the AF3 gate")
    C.announce(f"{cfg.target.name} af3_prescreen done",
               f"{len(passed)}/{len(rows)} champions pass AF3 "
               f"(iptm>={thr.iptm_min}, pae<={thr.chain_pae_max}, "
               f"plddt_b>={thr.plddt_min})", tags="checkered_flag")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
