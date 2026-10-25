#!/usr/bin/env python3
"""Round 1 (de novo): RFdiffusion -> MPNN/FastRelax -> AF2 -> QC -> funnel.

    python campaign/run_round1.py campaign/mouse_round1/round1.yaml --check
    python campaign/run_round1.py campaign/mouse_round1/round1.yaml --no-smoke

AF3 is not part of round 1; run campaign/validate_af3.py afterwards.
"""
from __future__ import annotations

import argparse
import os

import _campaign as C


def describe(cfg) -> None:
    from binder.rounds import round1
    t = cfg.target
    print(f"target        {t.structure}")
    print(f"host chain    {t.reference_chain}   copies {t.target_copies}")
    print(f"contig        [{cfg.binder.segment()} {round1._target_contig(cfg)}]")
    print(f"hotspots      {round1._hotspots_contig(t)}")
    print(f"recipe        {round1._recipe(cfg, cfg.binder.length)}")
    print(f"patch         {len(t.effective_patch())} residues")
    print(f"keep_assembly {t.keep_assembly}   species {t.species}")
    print(f"backbones     {cfg.round.backbones} x {cfg.round.mpnn_passes} passes "
          f"= {cfg.round.backbones * cfg.round.mpnn_passes} designs")
    print(f"mpnn          T={cfg.campaign.mpnn.temperature} "
          f"refine={cfg.campaign.mpnn.refine}")
    print(f"resources     gpus {cfg.resources.gpus}  "
          f"mpnn_workers {cfg.resources.mpnn_workers}  "
          f"af2_per_gpu {cfg.resources.af2_per_gpu}")
    print(f"results       {cfg.results_root}")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("config")
    ap.add_argument("--check", action="store_true",
                    help="resolve and verify inputs, print the plan, exit")
    ap.add_argument("--no-smoke", action="store_true")
    ap.add_argument("--lengths", nargs="*", type=int)
    ap.add_argument("--per-length", type=int)
    args = ap.parse_args(argv)

    cfg = C.load(args.config)
    C.require_inputs(cfg)
    describe(cfg)
    print("\nall inputs present")
    if args.check:
        return 0

    from binder.rounds import round1
    log = os.path.join(cfg.run_root, "round1", "run.log")
    notifier = C.setup(cfg, log=log)
    n_designs = cfg.round.backbones * cfg.round.mpnn_passes
    C.announce(f"{cfg.target.name} round1 start",
               f"{cfg.round.backbones} backbones (len {cfg.binder.length})\n"
               f"hotspots {cfg.target.hotspots}\n"
               f"{cfg.round.mpnn_passes}x [MPNN -> FastRelax -> MPNN], "
               f"AF2 on {n_designs} designs", tags="rocket")
    describe(cfg)

    # arm() sends the closing message whichever way this ends, with the tail of
    # the log on failure.
    with notifier.arm(f"{cfg.target.name} round1", log):
        result = round1.run(cfg, lengths=args.lengths,
                            per_length=args.per_length,
                            smoke=not args.no_smoke,
                            notify=lambda step: C.announce(
                                f"{cfg.target.name} {step}", f"{step} finished",
                                tags="white_check_mark"))
        stages = result.get("stages", {})
        failed = result.get("failed", [])
        n_pass = sum(len(s.get("pass", [])) for s in stages.values())
        n_short = sum(len(s.get("shortlist", [])) for s in stages.values())
        notifier.note(f"AF2 gate {n_pass}, shortlist {n_short}")
        notifier.note(f"base: {result.get('base')}")
        # run() carries on past a stage failure, so a dead stage must not read
        # as done: raise so arm() reports it as a failure.
        if failed:
            raise SystemExit(f"stage(s) {failed} failed; on the rest: AF2 gate "
                             f"{n_pass}, shortlist {n_short}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
