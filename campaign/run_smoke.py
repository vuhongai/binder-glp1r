#!/usr/bin/env python3
"""Smoke-only: run a small de novo batch in the smoke stage to validate a
config (length, hotspots, backbone filter) before committing to the full run.

    python campaign/run_smoke.py campaign/glp1r_round1/round1.yaml --num 16
    python campaign/run_smoke.py campaign/glp1r_round1/round1.yaml --check

Writes to outputs/campaign/<target>/round1/smoke_len<length>/ -- a separate
directory from the full stage1_len<length>, so it never collides with a full
run and can be re-run freely.  It applies the AF2 gate and funnel so the pae
distribution is visible at a glance: if the best pae_interaction is not
dropping toward the gate (10), the config is not ready for the full run.
"""
from __future__ import annotations

import argparse
import os

import _campaign as C


def describe(cfg) -> None:
    from binder.rounds import round1
    t = cfg.target
    print(f"target        {t.structure}")
    print(f"contig        [{cfg.binder.segment()} {round1._target_contig(cfg)}]")
    print(f"hotspots      {round1._hotspots_contig(t)}")
    print(f"binder len    {cfg.binder.length}")
    print(f"filter        {cfg.round.filter_backbones} "
          f"(extent<={cfg.round.filter_extent_max} "
          f"rg<={cfg.round.filter_rg_max} "
          f"hot_dist<={cfg.round.filter_hotspot_dist_max})")
    print(f"mpnn          T={cfg.campaign.mpnn.temperature} "
          f"refine={cfg.campaign.mpnn.refine}")
    print(f"resources     gpus {cfg.resources.gpus}  "
          f"af2_per_gpu {cfg.resources.af2_per_gpu}")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("config")
    ap.add_argument("--num", type=int, default=8,
                    help="backbones for the smoke stage (default 8)")
    ap.add_argument("--check", action="store_true",
                    help="resolve and verify inputs, print the plan, exit")
    args = ap.parse_args(argv)

    cfg = C.load(args.config)
    C.require_inputs(cfg)
    describe(cfg)
    print("\nall inputs present")
    if args.check:
        return 0

    from binder.rounds import round1
    length = cfg.binder.length
    base = os.path.join(cfg.run_root, "round1")
    root = os.path.join(base, f"smoke_len{length}")
    log = os.path.join(base, "run.log")
    notifier = C.setup(cfg, log=log)
    C.announce(f"{cfg.target.name} smoke start",
               f"{args.num} backbones (len {length})\n"
               f"hotspots {cfg.target.hotspots}\n"
               f"validation before the full run", tags="rocket")
    with notifier.arm(f"{cfg.target.name} smoke", log):
        res = round1.run_stage(cfg, root, f"stage1_len{length}", args.num,
                               length, round1._diffuse,
                               notify=lambda step: C.announce(
                                   f"{cfg.target.name} smoke {step}",
                                   f"{step} finished",
                                   tags="white_check_mark"))
        res["pass"] = round1.af2_pass(cfg, res["paths"])
        round1._structural_qc(cfg, res["paths"])
        res["shortlist"] = round1._funnel(cfg, res["paths"])
        n_pass = len(res.get("pass", []))
        n_short = len(res.get("shortlist", []))
        notifier.note(f"smoke AF2 gate {n_pass}, shortlist {n_short}")
        notifier.note(f"base: {base}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
