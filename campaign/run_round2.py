#!/usr/bin/env python3
"""Round 2 (basin refinement): partial diffusion around round-1 champions,
then the deep funnel (AF2 -> AF3 -> Rosetta -> composite ranking).

    python campaign/run_round2.py campaign/glp1r_round2/round2.yaml --check
    python campaign/run_round2.py campaign/glp1r_round2/round2.yaml
    python campaign/run_round2.py campaign/glp1r_round2/round2.yaml \
        --champions g0r0_677=<pdb> g1r0_36=<pdb>

Champions default to the round-1 shortlist (every AF2 gate-passer with an
extracted PDB). Pass --champions to pick a subset as tag=source_pdb pairs.

HOTSPOT REMAPPING (the round-2 fix).  A champion pose is chain A = binder
(1-90), chain B = target (renumbered 91-474).  RFdiffusion's hotspot list must
name the TARGET chain (B) with its renumbered positions.  Deriving it from the
config's reference_chain + hotspots (as the TNF round-2 did) silently points at
the binder's own residues and the partial diffusion runs unsteered.  This
driver therefore re-derives the hotspots for EVERY champion via
campaign/remap_hotspots.py -- mapping by alignment and verifying the residue
identity at each mapped position -- and refuses to start if the champions
disagree or any identity check fails.
"""
from __future__ import annotations

import argparse
import os

import _campaign as C
import remap_hotspots


def _apply_hotspot_remap(cfg, champions: list[tuple[str, str]]) -> str:
    """Compute the round-2 hotspots_contig (chain B, renumbered) for every
    champion, require them all to agree, and set it on the config.  Returns the
    contig string.  Raises SystemExit on any identity mismatch or disagreement.
    """
    target_pdb = cfg.target.structure
    hotspots = list(cfg.target.hotspots)
    if not hotspots:
        raise SystemExit("no target.hotspots set -- nothing to remap")

    seen: dict[str, str] = {}
    for tag, src in champions:
        contig, details = remap_hotspots.compute_hotspots_contig(
            src, target_pdb, hotspots)
        if contig in seen:
            if seen[contig] != tag:
                pass  # same contig from a different champion is fine
        else:
            seen[contig] = tag
        print(f"  {tag:16s} -> {contig}")
        for h, c, r, ident in details:
            print(f"      design {h:3d}  {c}{r} ({ident})")

    if len(seen) != 1:
        raise SystemExit(
            "the champions map to DIFFERENT hotspot contigs -- they are not all "
            "the same pose with the same offset, so a single hotspots_contig "
            "cannot steer them all:\n  " + "\n  ".join(sorted(seen)))
    contig = next(iter(seen))

    # The installed package must actually carry the field; if it silently drops
    # it the diffusion would run unsteered (the TNF round-2 bug).
    if not hasattr(cfg.target, "hotspots_contig"):
        raise SystemExit(
            "this binder package has no target.hotspots_contig field -- the "
            "remapped hotspots would be dropped and the diffusion unsteered")
    cfg.target.hotspots_contig = contig
    return contig


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("config")
    ap.add_argument("--champions", nargs="*",
                    help="tag=source_pdb pairs; default: the round-1 shortlist")
    ap.add_argument("--stage", help="round-1 stage to take champions from")
    ap.add_argument("--num", type=int, help="partial-diffusion backbones per champion")
    ap.add_argument("--check", action="store_true")
    args = ap.parse_args(argv)

    cfg = C.load(args.config)
    C.require_inputs(cfg)
    from binder.rounds import round2k

    specs = list(args.champions or [])
    if specs:
        champions = []
        for spec in specs:
            tag, _, src = spec.partition("=")
            if not src:
                raise SystemExit(f"--champions expects tag=source_pdb, got {spec}")
            src = src if os.path.isabs(src) else os.path.join(C.ROOT, src)
            if not os.path.isfile(src):
                raise SystemExit(f"champion PDB not found: {src}")
            champions.append((tag, src))
    else:
        champions = round2k.discover_champions(cfg, args.stage)
        if not champions:
            raise SystemExit(
                "no round-1 shortlist found -- round 1 must finish its AF2 gate "
                "first, or pass --champions tag=source_pdb")

    # Two champions on the same diffusion backbone spend a full round-2 budget
    # twice on nearly the same starting pose.
    def backbone_of_tag(tag: str) -> str:
        stem = tag.rsplit("_p", 1)
        return stem[0] if len(stem) == 2 and stem[1].isdigit() else tag

    seen_bb: dict[str, str] = {}
    clashing = []
    for tag, _ in champions:
        bb = backbone_of_tag(tag)
        if bb in seen_bb:
            clashing.append((seen_bb[bb], tag, bb))
        else:
            seen_bb[bb] = tag
    if clashing:
        print(f"WARNING: {len(clashing)} champion(s) share a backbone with "
              f"another; each costs a full "
              f"{(args.num or cfg.round.backbones) * cfg.round.mpnn_passes} "
              f"designs for nearly the same starting pose:")
        for first, dup, bb in clashing:
            print(f"  {bb}: {first} and {dup}")
        print(f"  {len(champions)} champions over {len(seen_bb)} backbones. "
              f"Pass --champions explicitly (one design per backbone).")
        raise SystemExit(1)

    # THE FIX: remap the hotspots onto the champion pose before anything runs.
    num = args.num or cfg.round.backbones
    print(f"champions     {len(champions)} over {len(seen_bb)} backbone(s)")
    print("hotspot remap (chain B, renumbered):")
    contig = _apply_hotspot_remap(cfg, champions)
    print(f"hotspots_contig set to: {contig}")

    print(f"per champion  {num} partial-diffusion backbones x "
          f"{cfg.round.mpnn_passes} MPNN passes = {num * cfg.round.mpnn_passes} designs")
    print(f"partial_t     {cfg.campaign.partial_t}")
    print(f"funnel        AF2 gate -> AF3 ({cfg.target.species}) -> Rosetta -> rank")
    print(f"results       {os.path.join(cfg.run_root, 'round2')}")
    if args.check:
        return 0

    log = os.path.join(cfg.run_root, "round2", "run.log")
    notifier = C.setup(cfg, log=log)
    C.announce(f"{cfg.target.name} round2 start",
               f"{len(champions)} champions x {num} backbones x "
               f"{cfg.round.mpnn_passes} passes\nhotspots {contig}",
               tags="rocket")
    with notifier.arm(f"{cfg.target.name} round2", log):
        result = round2k.run(cfg, champions, num)
        failed = result.get("failed", [])
        notifier.note(f"{len(result.get('champions', {}))} champions done, "
                      f"{len(failed)} failed")
        if failed:
            raise SystemExit(f"champions failed: {failed}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
