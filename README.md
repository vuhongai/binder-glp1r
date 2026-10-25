# GLP1R agonist binder — 2-round design campaign

De novo binders for human GLP1R that both bind and activate (Gs/cAMP).
Two rounds, mirroring `/DATA/binder-hmTNF`:

## Status (2026-10-19)

- **Round 1 (90-mer) gave 0 leads** (3056 designs, best pae_interaction 11.5,
  gate 10). Diagnosis: the 90-mers (50–73 Å) are 3.6–4× the 18 Å hotspot patch,
  so they drape over the whole TM domain (108–132 residues contacted) instead of
  a local interface → AF2 can't be confident. See `memory.md`.
- **Rerunning as 60-mers** (25–40 Å, wraps the patch locally), all 6 hotspots,
  1000 backbones × 8 = 8000 designs. Deadline **2026-10-25**.
- Workflow: `campaign/run_smoke.py --num 16` (validate pae drops) →
  `run_round1.py --no-smoke` (full) → `run_round2.py` on the champions.

- **Round 1 (de novo):** RFdiffusion → MPNN/FastRelax → AF2 → gate → funnel.
- **Round 2 (basin refinement):** partial diffusion around round-1 champions →
  MPNN → AF2 → gate → AF3 → Rosetta → composite ranking.

Engine: the `/DATA/binder` package. This folder only holds the target, the
configs, and the thin drivers.

## Target and hotspots

- **Structure:** 6X18 receptor, author chain R (the active-state, GLP-1-bound
  human receptor — the organizer's structural reference). The **ECD is trimmed**
  (it sits 34.7 Å from the hotspots and competes for the binder), so the
  diffusion target is the **TM domain + C-tail, design 73–384** (UniProt
  101–423) in `data/targets/glp1r/receptor_tmdomain.pdb`. The full receptor
  (design 1–384) is kept at `receptor.pdb` for reference.
- **Numbering:** design number = residue number in the PDB. The design →
  UniProt map is **not a constant offset** (4 loops are unmodelled); use
  `data/targets/glp1r/design2up_tmdomain.json` as the authoritative map. Full
  table: `data/targets/glp1r/target_report.csv`.
- **Hotspots** (design / UniProt / residue), all in the TM domain:

  | design | UniProt | residue | region (P1/P2 nomenclature) |
  |---|---|---|---|
  | 264 | 298 | T | ECL2 |
  | 265 | 299 | R | ECL2 |
  | 266 | 300 | N | ECL2 |
  | 333 | 372 | D | ECL3/TM7 |
  | 341 | 380 | R | ECL3/TM7 |
  | 345 | 384 | L | ECL3/TM7 |

  The context.md "primary seed pool" (T298/R299/N300) plus the "complementary
  activation-facing candidates" (R380, with D372 and L384 as support). **All six
  directly contact GLP-1 in 6X18** (verified against P1 Table S2; see
  `glp1_contacts.json`). The ECL2 trio is compact (5–9 Å); D372/R380/L384 sit
  ~16 Å away in the ECL3/TM7 region.

### Known deviations from wild-type (do not treat as hotspots)
- **L260F:** the 6X18 construct has Phe at UniProt 260 (Leu in wild-type), in
  ECL1. It is the single deviation in the 383/384 Needleman-Wunsch match to
  UniProt P43220. It is not a hotspot and is not contacted by GLP-1.
- **Construct uncertainty:** the organizer's assay construct (13944-H18H-DD,
  UniProt 24–432) has undisclosed thermostabilizing mutations. This design uses
  the wild-type 6X18 sequence (the organizer's structural reference). The
  trimmed target keeps the C-tail to UniProt 423 (design 305–384), beyond the
  organizer's 24–432 reference; it is on the intracellular face, far from the
  hotspots, and is kept as diffusion context.
- **ECD trimmed:** the N-terminal ECD (design 1–72, UniProt 29–100) is removed
  from the diffusion target. It is 34.7 Å from the hotspots and its large
  surface was drawing the de novo backbones off the TM face (the smoke test
  gave elongated backbones, pae_interaction 13–28). AF2/AF3 therefore validate
  binder+TM (they do not see the ECD); since the binder binds the TM face ~10 Å
  from the hotspots, an ECD clash is unlikely but is a small blind spot.

## The round-2 hotspot fix (important)

A round-2 champion pose is **chain A = binder (1–90), chain B = target
(renumbered 91–474)**. RFdiffusion's hotspot list must name the **target chain
(B)** with its renumbered positions. Deriving it from the config's
`reference_chain` + `hotspots` (what the TNF round-2 did) silently points at
the **binder's own residues**, and the partial diffusion runs unsteered.

`campaign/remap_hotspots.py` fixes this: it aligns each champion's target chain
to the original target, maps every hotspot through the alignment (never a fixed
offset), and **verifies the residue identity at each mapped position**.
`run_round2.py` calls it for every champion, requires them all to agree, and
refuses to start on any identity mismatch or disagreement. For this target the
remapped list is `[B354,B355,B356,B423,B431,B435]` (each design number + 90).

Self-test (CPU only, no GPU):

    python campaign/remap_hotspots.py --self-test \
        data/targets/glp1r/receptor_tmdomain.pdb 264,265,266,333,341,345

## The pre-MPNN backbone filter

`round.filter_backbones` (in the binder package, off by default) measures each
backbone's CA extent, rg and min binder→hotspot distance and moves the failures
to `<stage>/diffusion/rejected/` before MPNN/FastRelax/AF2 spend time on them.

The de novo backbones are consistently **elongated** (extent ~64–73 Å, rg ~20 Å;
a compact 90-mer is ~30–40 Å) yet sit on the hotspot region (hot_dist ~4–10 Å).
So the **docking mode** (hot_dist) is the active gate and the shape limits are
kept loose (extent ≤ 80 Å, rg ≤ 25 Å, hot_dist ≤ 12 Å — tune in `base.yaml`);
a hard shape gate would reject the whole batch. A **floor** of `max(4, n//10)`
promotes the most compact failures back in, so the filter can never empty the
stage (the first re-run crashed on exactly that: 0/8 kept → empty `bb.silent`).

## How to run

The LLM server holds both GPUs; **stop it before launching** (the drivers hold
for a free card and notify). All commands run from the repo root
(`/DATA/binder-glp1r`) with the `binder` env.

### Round 1

    PY=/home/avuhong/anaconda3/envs/binder/bin/python
    # 0. validate a config change cheaply (separate smoke_len60 dir, AF2 gate
    #    applied so the pae distribution is visible before the full run):
    $PY campaign/run_smoke.py campaign/glp1r_round1/round1.yaml --num 16
    # 1. verify inputs + plan (CPU only):
    $PY campaign/run_round1.py campaign/glp1r_round1/round1.yaml --check
    # 2. launch the full run (smoke already done, so skip it):
    $PY campaign/run_round1.py campaign/glp1r_round1/round1.yaml --no-smoke

Default scale: **1000 backbones × 8 MPNN passes = 8000 designs**. Scale with
`round.backbones` in `round1.yaml` or `--per-length <n>`.

### Round 2 (after round 1 has a shortlist)

    $PY campaign/run_round2.py campaign/glp1r_round2/round2.yaml --check
    $PY campaign/run_round2.py campaign/glp1r_round2/round2.yaml

Champions default to the round-1 shortlist. Pick a subset with
`--champions tag=source_pdb ...` (e.g. only the AF3-confirmed, clash-free
designs). Scale with `round.backbones` (16/champion) or `--num <n>`.

Smoke-test the round-2 driver (remap + plan, no GPU) before round 1 finishes,
using the bundled synthetic champion (90-res binder + target renumbered +90):

    $PY campaign/run_round2.py campaign/glp1r_round2/round2.yaml --check \
        --champions test0=data/test/synthetic_champion.pdb

### Submission (any time, CPU only)

`campaign/submit.py` builds the Proteinbase upload from **whatever has finished** —
safe to run on a half-done round 2, and safe to re-run (it caches the structural
audit and never overwrites `methodology.md`).

    $PY campaign/submit.py --n 40          # Track 1; use --n 20 for Tracks 2/3

Writes to `outputs/submission/`: `submission.csv` (name, sequence,
molecule_class), `submission_annotated.csv` (every metric used),
`structure_audit.csv` (per-design hotspot contacts + epitope size) and
`selection_report.md` (what was pooled, filtered and why). Champions still mid-run
are reported as skipped rather than silently dropped.

Selection is **agonism-directed**: ECL3/TM7 engagement (D372/R380/L384) is a hard
requirement wherever an AF3 model is reachable, since a binder that does not
activate scores zero. Run it on the machine holding the AF3 outputs — `af3/` is
gitignored, and without those models the hotspot terms fall back to neutral.

## Key parameters (base.yaml)

- binder length **90**; MPNN **temperature 0.15**, `refine: false`
  (round-1.5 protocol — no argmax refine, the source of the Ala bias).
- **backbone filter on:** extent ≤ 80 Å, rg ≤ 25 Å, hot_dist ≤ 12 Å, floor
  `max(4, n//10)`.
- partial_t **10**; max_winners **6**.
- GPUs `[0, 1]`; mpnn_workers 64; af2_per_gpu 2; rosetta_parallel 16.
- ntfy topic `design_glp1r_binder`.

## Layout

    data/targets/glp1r/   receptor_tmdomain.pdb (+.contig, design 73-384) =
                          the diffusion target; receptor.pdb = full 1-384
                          (reference); design2up_tmdomain.json,
                          glp1_contacts.json, target_report.csv
    data/fasta/           GLP1R_HUMAN_P43220.fasta (AF3 target MSA)
    data/test/            synthetic_champion.pdb (round-2 smoke test)
    campaign/
      base.yaml           shared config (target/binder/mpnn/thresholds/resources)
      glp1r_round1/       round1.yaml (extends base)
      glp1r_round2/       round2.yaml (extends base)
      _campaign.py        glue: config load, path resolution, notify, stop
      run_round1.py       round-1 driver
      run_round2.py       round-2 driver (wires in the hotspot remap)
      remap_hotspots.py   round-2 hotspot remapping + self-test
    outputs/campaign/     results (created at run time)
