# GLP1R agonist binder — memory

Code versions, bugs, campaign parameters, candidate identifiers, outcomes.
Scientific rationale lives in `context.md`; run instructions in `README.md`.

## Campaign (built 2026-10-17)

- **Goal:** de novo 90-mer human GLP1R agonist (bind + Gs/cAMP activation).
- **Engine:** `/DATA/binder` package; 2 rounds (de novo → basin refinement),
  mirroring `/DATA/binder-hmTNF`.
- **Target:** 6X18 receptor author chain R, **ECD trimmed** → TM domain + C-tail,
  **design 73–384** (`data/targets/glp1r/receptor_tmdomain.pdb`, contig
  `A73-384`). Full receptor (design 1–384) kept at `receptor.pdb` for reference.
  Design→UniProt is **not a constant offset** (4 unmodeled loops); authoritative
  map = `data/targets/glp1r/design2up_tmdomain.json`.
- **Hotspots (design/UniProt/residue):** 264/T298, 265/R299, 266/N300,
  333/D372, 341/R380, 345/L384. = context.md primary seed (T298/R299/N300) +
  complementary (R380 + D372/L384 support). All six contact GLP-1 in 6X18
  (validated vs P1 Table S2). All in the TM domain.
- **Params:** binder **60 aa** (was 90 — see "Round 1 (90-mer) outcome"); MPNN
  T=0.15, refine=false (no argmax refine → no Ala bias); **backbone filter on**
  (60-mer: extent≤50, rg≤18, hot_dist≤12); partial_t 10; max_winners 6; GPUs
  [0,1]; ntfy `design_glp1r_binder`.
- **Scale:** round 1 = 1000 backbones × 8 passes = 8000 designs; round 2 =
  16 backbones/champion × 8 passes.

## Mapping validation (done, CPU)

- 6X18 receptor = auth chain R (label F, entity 6). UniProt P43220 FASTA is
  **463 aa** (not 438). Receptor proper = entity seq_id 33–472 = UniProt
  24–463; `UniProt = seq_id − 9`.
- Modeled receptor = 384 residues; Needleman-Wunsch to UniProt = 383/384 match.
  Single deviation = **L260F** (modeled Phe vs UniProt 260 Leu, ECL1) — a
  construct mutation, not a hotspot, not GLP-1-contacted.
- GLP-1 contact map (receptor within 4 Å of GLP-1 auth chain P) reproduces P1
  Table S2 almost exactly → validates mapping + template.
- Apo SASA (Å²): T298 40.9, R299 87.2, N300 60.8, R380 78.6, L384 67.6,
  D372 42.0, W33 142.5, E68 122.8, R121 91.3, W91 83.8.
- 6LN2 (inactive, auth chain A) alignment unreliable (rubredoxin insertion +
  11 mutations → frame shift, 142 Å superposition) → recorded as limitation,
  not used for the hotspot decision.

## Round-2 hotspot bug — fixed

- **Bug:** champion pose is chain A = binder (1–90), chain B = target
  (renumbered 91–474). Deriving `hotspots_contig` from `reference_chain`+
  `hotspots` (TNF round 2) points at the **binder's own** residues → unsteered
  diffusion.
- **Fix:** `campaign/remap_hotspots.py` aligns each champion's target chain to
  the original target, maps hotspots through the alignment (no fixed offset),
  and verifies residue identity at each mapped position. `run_round2.py` runs
  it for every champion, requires agreement, refuses on mismatch/disagreement.
- **Remapped list for this target (60-mer champions, verified 2026-10-23):**
  `[B252,B253,B254,B321,B329,B333]` (design + 88 offset). All 42 round-1
  survivors map to this one contig with verified residue identity.
- Self-test passes (happy path + wrong-identity refusal + disagreeing-champion
  refusal). `run_round1.py --check` and `run_round2.py --check` both pass.

## Smoke test (2026-10-17) → ECD trim + backbone filter

- **Symptom:** 8 de novo backbones all elongated (CA extent 57–107 Å, rg
  15–33 Å; compact 90-mer ~30–40 / ~12–15 Å) and binder centroid 13–32 Å from
  the hotspots; AF2 pae_interaction 13–28 (a binder is <10), target_aligned_rmsd
  up to 69 Å.
- **2-domain analysis (CA centroids, trimmed-receptor frame):** ECD (design
  1–72, UP 29–100) y≈177; TM (73–304, UP 101–343) y≈134; C-tail (305–384,
  UP 344–423) y≈117. ECD↔TM 43 Å apart. **All hotspots are in the TM domain**
  (10.8 Å from TM centroid); ECD is 34.7 Å from the hotspots. Whole-receptor
  extent 123 Å.
- **Decision (user-confirmed):** (1) trim the ECD → diffusion target is TM +
  C-tail (design 73–384); (2) add D372 + L384 hotspots (6 total); (3) add a
  pre-MPNN backbone filter.
- **Backbone filter (new, in `/DATA/binder`):** `round.filter_backbones`
  (off by default) measures CA extent, rg, min binder→hotspot distance; moves
  failures to `<stage>/diffusion/rejected/`. Called in `round1.run_stage`
  between diffusion and packing. Enabled in `base.yaml`.
- **Filter crash fix (2026-10-17):** first re-run (trimmed target) the filter
  rejected all 8 smoke backbones (extent 64–73 Å > 50 limit) → 0 backbones →
  empty `bb.silent` → `silentfrompdbs` rc=1 crash. Two fixes: (1) loosened the
  shape limits (extent≤80, rg≤25; hot_dist≤12) because the backbones are
  consistently elongated (64–73 Å, mean 70) yet on-target (hot_dist 3.5–9.7 Å)
  — the docking mode is the real gate, not shape; (2) added a **floor**
  `max(4, n//10)` that promotes the most compact failures back in, so the
  filter can never empty the stage. Stale smoke outputs deleted for a clean
  re-run.

## gemmi 0.7.5 / Python 3.14 notes (non-standard API)

- No `.models`/`.chains`; use `st[0]` (model), iterate model→chains.
- `residue.label_seq` (settable), `residue.seqid` (SeqId, what write_pdb uses),
  `residue.name`, `residue.get_ca()`, iterate residue→atoms, `atom.occ`,
  `atom.b_iso`, `atom.element.name`.
- `SeqId(int,str)` constructor broken → **write PDB manually** to renumber.
- The receptor PDB uses a non-standard column layout: resName cols 17–19,
  chain col 21, resSeq cols 22–25 (the binder package parses the same).

## Round 1 (90-mer) outcome + diagnosis (2026-10-19)

- **Result: 0 leads.** 382 backbones × 8 passes = 3056 designs; AF2 gate
  (pae<10, rmsd<2, plddt_b>80) **0/3056 pass**. Best pae_interaction **11.52**
  (2 designs <12, 100 <15, median 27). binder_aligned_rmsd median 1.8 and
  plddt_binder median 75 are fine — the *fold* is not the problem.
- **Diagnosis (geometry, top-20 poses measured vs the hotspot patch):** the
  binders bind the *right place* (min binder→target 0.2–1.0 Å; all 6 hotspots
  within 1.3–3.5 Å) but they are **too big**. Every top-20 90-mer contacts
  **108–132 target residues spanning design 73→355** — they drape over the
  *entire* TM domain, not a local interface. The 90-mer backbones are
  **50–73 Å** in extent vs the **18 Å** hotspot patch (ECL2 trio + ECL3 cluster,
  14 Å apart) → 3.6–4× the patch → flat, non-specific interface → AF2 reports
  low interaction confidence (pae 11–13) and drags the target
  (target_aligned_rmsd median 33).
- **Compounding:** only 3056 designs (TNF used 16000); and TNF is a *trimer*
  (deep pocket) while GLP1R is a *monomer* with a flat face (harder). TNF got
  20 hits (0.125%) at 16000; GLP1R at 3056 with an over-long binder got 0.
- **Decision (user-confirmed 2026-10-19):** rerun as **60-mers** (a 60-mer is
  25–40 Å, only 1.4–2.2× the patch → wraps it locally → tighter, confident
  interface). Keep **all 6 hotspots** (ECL3/R380 needed for agonism). Scale up
  to **1000 backbones × 8 = 8000 designs** (60-mers are cheaper). Tightened
  backbone filter for 60-mers: extent ≤50, rg ≤18 (in base.yaml).
- **Deadline pressure:** challenge deadline **2026-10-25** (~6 days). Plan:
  60-mer smoke (16 bb) → validate pae drops → full 1000-bb round 1 (~1.5–2 d)
  → round 2 on champions → submit.
- **New driver:** `campaign/run_smoke.py` — smoke-only (separate `smoke_len60`
  dir, never collides with the full `stage1_len60`), applies the AF2 gate so
  the pae distribution is visible before committing to the full run. The engine
  has no smoke-only mode (run() always does smoke-then-full; and the resume
  logic re-diffuses the full `num` when `have < num`, so a small run cannot be
  extended in-place).

## Outcomes

- Round 1 (90-mer, 3056 designs): **0 AF2-gate passes, 0 leads.** Best
  pae_interaction 11.52. Root cause: 90-mers too long → drape over the whole
  TM domain (108–132 residues) instead of a local 18 Å patch interface.
- Round 1 (60-mer smoke, 16 bb, 2026-10-19): 6/16 backbones survived the
  filter (extent≤50 rejected 60%, incl. 3 at 51–53 Å that are essentially
  compact). 48 designs; best pae_interaction **10.42** (vs 11.52 for the 90-mer),
  target_aligned_rmsd median **9** (vs 33) — the 60-mer is a real improvement
  (more local interface, stable target) but marginal (best pae just above the
  gate of 10, 0 pass on 48 designs).
- **Filter loosened** (base.yaml) to extent≤60, rg≤20: keeps the near-compact
  51–60 Å backbones (the good ones) while still rejecting the truly elongated
  67–86 Å ones. ~56% of backbones now survive.
- Round 1 (60-mer full, 1000 bb × 8 = 8000 designs): launched 2026-10-19
  21:32 in tmux `glp1r_full`, holding for a free GPU; LLM killed to start it.
  Will ntfy on completion (~1.5–2 d).

## Round 1 (60-mer full) outcome — DONE 2026-10-23 05:22

- **Result: 91 AF2-gate passes (1.1%), 42 funnel survivors.** 634 backbones
  survived the filter (extent≤60/rg≤20/hot_dist≤12); 8000 designs scored.
  Gate (pae<10, binder_aligned_rmsd<2, plddt_binder>80): **91/8080 pass**.
  Best pae_interaction **8.07** (`g0r0_451_p7`, rmsd 0.39, plddt_b 86.7).
  Funnel (pae<10, one per backbone): **42 survivors**, pae 8.07–9.85, 42
  distinct backbones. Big jump vs 90-mer (0/3056) and smoke (0/48).
- **Stop/restart validity (audited 2026-10-23): DATA IS VALID.** The run was
  killed/restarted 5× (10-19 21:32, 10-20 14:35, 10-20 21:07, 10-21 09:15,
  10-21 15:00, 10-21 21:52→done). `cautious=True` re-diffused only missing
  indices; run #3's filter rejected 618 backbones (moved to `diffusion/
  rejected/`), so later restarts re-diffused those 618 indices → 252 of the
  regenerated passed the final filter. Final kept set = 634 = 382 original +
  252 regenerated (all legitimate new de novo backbones). Stale copies in
  `diffusion/rejected/` are harmless (pipeline never reads that dir).
- **Stale MPNN parts (harmless, noted):** block dirs b1/b2/b3 still hold 64
  old part files (254 backbones) from run #3 that the final AF2 also scored
  (same tags, older MPNN seqs) → duplicate AF2 samples; scores.csv = 8080 rows
  = 5072 current (634×8) + 3008 stale-part samples. All 91 shortlist PDBs
  verified to match their tag's AF2 structure; 13 passes are on run #3 MPNN
  seqs but the PDB is self-contained (valid for round 2). Do NOT re-run
  `run_round1.py` on this stage (`.sc` files append, stale parts re-scored).

## Round 2 prep (verified 2026-10-23, CPU)

- **Champions:** the 42 funnel survivors = exactly one design per backbone
  (the best-pae passer of each). `run_round2.py` refuses auto-discovery here
  (91 gate-passers span only 42 backbones → multiple designs/backbone), so
  pass `--champions` explicitly. List saved to
  `campaign/champions_round2.txt` (42 lines, `short_id=<rel pdb path>`).
- **`run_round2.py --check` with the 42 champions PASSES:** all 42 remap to the
  single contig `[B252,B253,B254,B321,B329,B333]` with verified identity.
  Scale: 42 × 16 partial-diffusion backbones × 8 MPNN passes = 5376 designs;
  partial_t 10; funnel AF2 gate → AF3 (human) → Rosetta → rank; results under
  `outputs/campaign/glp1r/round2`.
- **Launch (blocked on user):** needs a free GPU — both A30s are held by the
  LLM server (~24 GB each); the driver waits for a free card. Deadline
  2026-10-25 23:59 AoE (~2 d). 42 champions ≈ 2 d wall (tight); a top-15
  subset ≈ 1–1.5 d leaves buffer for AF3/Rosetta + submission.
- **Champions are NOT AF3-validated** — round-1 funnel used only the AF2
  initial-guess gate (pae<10, rmsd<2, plddt_b>80). AF3 is the blind check that
  round 2 applies to the new designs. `campaign/af3_prescreen.py` runs blind
  AF3 on the champion binder sequences (reuses the round-2 MSA cache) to narrow
  the set; needs a free GPU (~6.5 h). Top-15 list saved to
  `campaign/champions_round2_top15.txt`.

## AF3 pre-screen of the 42 champions (done 2026-10-23) + lead selection

- **Strict AF3 gate (iptm≥0.80, pae_TB≤10, plddt_b≥80): 0/42 pass.** The
  plddt_b≥80 term is the killer and is too tight for de novo binders (no MSA →
  AF3 under-confident on the binder alone; every champion's plddt_b is 62–79).
  The integrated **iptm** is the right readout.
- **iptm≥0.80 (the standard high-confidence-complex cut): 12/42 leads.** All
  also have pae_TB 1.67–2.27 (very confident interface) and plddt_b 72.6–79.1.
  Selected as the round-2 champions → `campaign/champions_round2_af3.txt`.
  `run_round2.py --check` passes: 12 over 12 backbones, all remap to the single
  contig `[B252,B253,B254,B321,B329,B333]`. 12×16×8 = 1536 designs.
- **Key finding: AF3 iptm and AF2 pae are poorly correlated here** — only 4 of
  the 12 AF3 leads are in the AF2-pae top-15. The AF3 pre-screen (the user's
  instinct) selects a different, more confident set than the AF2 ranking.
- The 12 leads (short_id, iptm, pae_TB, plddt_b): g1r0_399_p1 (0.84/1.67/76.2),
  g0r0_199_p5 (0.84/1.78/72.6), g0r0_299_p7 (0.83/1.67/74.8), g1r0_421_p4
  (0.83/1.78/77.9), g1r0_298_p5 (0.83/1.89/79.1), g1r0_269_p3 (0.82/2.01/74.3),
  g0r0_142_p6 (0.82/2.06/73.0), g1r0_248_p6 (0.82/2.08/72.9), g0r0_27_p3
  (0.82/2.14/72.7), g0r0_356_p2 (0.82/2.24/77.4), g1r0_366_p3 (0.81/2.19/72.6),
  g1r0_120_p7 (0.80/2.27/73.0).

## Resume hardening (2026-10-23, in `/DATA/binder`)

User will stop/resume the round-2 run; made every stage resumable safely:
- **`af2.merge_scores` dedups by design tag** (keeps the most recent
  prediction): a resumed run re-appends to the Rosetta `.sc` score files, so
  without this a tag appears twice (old + new sample) and the gate/funnel judge
  a stale row.
- **`round1.run_stage` writes `af2_stage_done`** (after the score merge) and
  **`round2k.run_champion` skips diffusion/MPNN/AF2 when it exists**: a resumed
  multi-champion run re-enters every champion from the top; without this a
  finished champion re-runs its whole FastRelax/AF2 pipeline (hours) before the
  cheap resumable gate/AF3/Rosetta could skip it.
- **`round1.run_stage` removes stale `part_*.silent` (+ their `*.mpnn.silent`)**
  from the block dirs that aren't in the current deal: a changed backbone set
  (filter rejected → re-diffused) otherwise leaves old chunks that get MPNN'd
  and AF2'd as duplicate designs (the round-1 stale-part hazard).
- Already resumable (no change): partial diffusion (cautious skip + count/
  recipe), AF2 (per-chunk `-checkpoint_name`), AF3 (per-design marker + MSA
  cache), Rosetta (per-pose JSON skip), gate/QC/funnel/finalize (idempotent).
- Resume = re-run the SAME command. `run_round1.py --check`, `run_round2.py
  --check`, and `af3_prescreen.py --check` all pass after the edits.
