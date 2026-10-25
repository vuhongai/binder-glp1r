# De novo GLP-1R agonist binders — methodology

60-residue de novo mini-binders designed against the active-state human GLP-1R
transmembrane face, targeting the peptide-agonist contact surface so that binding
is coupled to receptor activation rather than to inert surface engagement.

---

## 1. Target definition and hotspot rationale

**Structural reference.** PDB **6X18**, author chain R — the active-state,
GLP-1-bound human receptor. Choosing the agonist-bound state is deliberate: a
binder selected against the inactive conformation would at best be an antagonist.

**Mapping was validated before any design was run**, because every downstream
hotspot number depends on it:

- Needleman–Wunsch alignment of the modelled receptor (384 residues) to UniProt
  P43220 gives **383/384 identity**. The single deviation is **L260F** in ECL1, a
  construct mutation — not a hotspot and not contacted by GLP-1.
- The receptor↔GLP-1 contact map computed from 6X18 (4 Å heavy-atom) reproduces
  the published contact table almost exactly, which validates both the mapping
  and the template.
- The design→UniProt correspondence is **not a constant offset** (four loops are
  unmodelled), so an explicit per-residue map is used throughout rather than
  arithmetic.

**Hotspots — six residues, all of which directly contact GLP-1 in 6X18:**

| UniProt | Residue | Region | Apo SASA (Å²) |
|---|---|---|---|
| 298 | Thr | ECL2 | 40.9 |
| 299 | Arg | ECL2 | 87.2 |
| 300 | Asn | ECL2 | 60.8 |
| 372 | Asp | ECL3/TM7 | 42.0 |
| 380 | Arg | ECL3/TM7 | 78.6 |
| 384 | Leu | ECL3/TM7 | 67.6 |

The rationale is mechanistic, not geometric. The **ECL2 trio (T298/R299/N300)** is
the compact core of the GLP-1 contact surface (5–9 Å spread). The **ECL3/TM7
cluster (D372/R380/L384)** sits ~16 Å away and is the activation-facing side: the
conformational coupling that drives Gs engagement runs through the ECL3/TM7
region, so a binder that buries only ECL2 is a plausible *binder* but an
implausible *agonist*. Both sub-clusters are therefore required, and together
they define an **~18 Å epitope patch** — a number that went on to govern binder
length (§2).

**The extracellular domain (ECD) is excluded from the diffusion target.** Its
centroid lies **34.7 Å** from the hotspots, and in a pilot run its large surface
drew backbones away from the TM face entirely. The design target is therefore the
**TM domain + C-tail (UniProt 101–423)**; the full receptor is retained for
reference.

---

## 2. Binder design

Two rounds: de novo generation, then basin refinement around validated leads.
RFdiffusion → ProteinMPNN (T = 0.15, no argmax refine) + Rosetta FastRelax →
AlphaFold2 initial-guess → gate.

**Round 1 at 90 aa failed, and the failure was informative.** 3,056 designs,
**0 passed** the AF2 gate (best pae_interaction 11.5). Diagnosis by direct
measurement of the top poses: the binders docked onto the *right place* (all six
hotspots within 1.3–3.5 Å) but each contacted **108–132 target residues spanning
the whole TM domain**. At 50–73 Å of CA extent against an 18 Å patch, a 90-mer
cannot form a local interface — it drapes, and AF2 correctly reports low
confidence in a flat, non-specific contact.

**Round 1 was re-run at 60 aa** (25–40 Å ≈ 1.4–2.2× the patch, so the binder wraps
it locally). A pre-MPNN backbone filter (CA extent ≤ 60 Å, Rg ≤ 20 Å, binder→hotspot
≤ 12 Å, with a floor that prevents the stage from emptying) removed elongated
backbones before expensive stages. **8,000 designs → 91 AF2-gate passes (1.1%)**
(pae_interaction < 10, binder_aligned_rmsd < 2 Å, pLDDT_binder > 80), best
pae_interaction **8.07**, reducing to **42 survivors on 42 distinct backbones**.

**Blind AF3 re-prediction of the 42** (no initial guess, target MSA only) was used
as an orthogonal check. Two findings shaped the rest of the campaign:

1. A strict gate including **pLDDT_binder ≥ 80 passed 0/42** and was discarded as
   inappropriate: a de novo binder has no MSA, so AF3 is systematically
   under-confident about the binder chain in isolation. The **integrated ipTM** is
   the meaningful interface readout, and **12 designs reached ipTM ≥ 0.80**.
2. **AF3 ipTM and AF2 pae_interaction are poorly correlated on this target** —
   only 4 of the 12 AF3 leads fall in the AF2 pae top-15. Ranking on either metric
   alone would have discarded half of the confident complexes.

**Round 2** applies partial diffusion (partial_t = 10) around those 12 leads,
6/16 backbones × 8 MPNN sequences each, then AF2 gate → AF3 → Rosetta
InterfaceAnalyzer.

---

## 3. Binder selection

Selection is **agonism-directed**, because a binder that does not activate scores
nothing under the challenge rules. Implemented in `campaign/submit.py`.

**Hard requirements:** no AF3 inter-chain clash; AF3 ipTM ≥ 0.80; unique sequence;
length within 10–250 aa.

**Hotspot engagement, measured per design on the AF3 complex** (heavy-atom
contacts ≤ 4.5 Å, target chain read in UniProt numbering):

- **Full 6/6 hotspot coverage is required** — every submitted design contacts all
  three ECL2 residues *and* all three ECL3/TM7 residues. This is the §1
  activation-facing requirement enforced as a filter rather than an aspiration,
  and it is a real discriminator: only 42 of 167 otherwise-eligible designs
  achieve it.
- It is also a genuine trade-off, taken deliberately. Hotspot coverage and ipTM
  rank *differently* on this target — the highest-ipTM designs are not the ones
  with the most complete hotspot engagement — so requiring 6/6 costs roughly
  0.02–0.03 ipTM at the top of the list. Under rules where a non-activating
  binder scores zero, mechanism was given priority over interface confidence.

**Interface locality.** Epitope size is measured for every design as a direct
guard against the round-1 failure mode. For these 60-mers the measured epitope is
**median 47 target residues (range 36–73)**, versus **108–132** for the failed
90-mers — confirming the length change achieved a local interface. Designs with
inflated epitopes are penalised and those above 85 residues rejected.

**Physics-based checks** (Rosetta InterfaceAnalyzer): shape complementarity,
interface hydrogen bonds, packstat, and ddG. ddG is used **normalised per 1,000 Å²
of buried interface**, not raw — on a flat, oversized interface a large raw ddG is
a symptom of the drape pathology rather than evidence of a better binder.

**Pose convergence is a hard requirement.** The binder RMSD between the AF2 and
AF3 models, after superposing on the target, must be **below 4 Å**. Across the
eligible pool this metric has median 4.2 Å and reaches 10.4 Å, so a design can
carry high ipTM and full hotspot coverage while the two methods disagree about
where the binder actually sits; those were excluded rather than ranked down.

**Diversity.** Picks are capped per parent lead and per sequence-identity cluster
(< 70% pairwise). With an expected wet-lab hit rate of a few percent, a submission
of near-identical sequences is one experiment repeated, not N experiments; the set
deliberately spans distinct backbones and binding modes.

---

<!-- BEGIN AUTOGENERATED LEAD TABLE -- rewritten by campaign/submit.py -->

## 5. The 20 submitted designs

Ranked by the composite described in §3. Every metric produced by the pipeline is listed, grouped by its source.

### 5.1 Interface confidence and pose (AlphaFold2 / AlphaFold3)

AF2 is the initial-guess prediction used as the design gate; AF3 is a blind re-prediction (target MSA only, no initial guess).

| # | design | parent lead | AF2 pae_int | AF2 rmsd_b | AF2 pLDDT_b | AF3 ipTM | AF3 PAE_int | AF3 PAE_b | AF3 pLDDT_b | AF2/AF3 rmsd | clash |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | g1r0_1_p4 | g0r0_199_p5 | 7.29 | 0.48 | 87.5 | 0.86 | 1.75 | 13.9 | 76.5 | 3.60 | 0 |
| 2 | g1r0_5_p4 | g0r0_299_p7 | 8.69 | 0.56 | 85.6 | 0.86 | 1.58 | 12.3 | 79.9 | 2.59 | 0 |
| 3 | g1r0_1_p3 | g1r0_399_p1 | 8.47 | 0.54 | 85.9 | 0.85 | 1.48 | 12.5 | 78.0 | 2.73 | 0 |
| 4 | g0r0_5_p7 | g0r0_299_p7 | 9.53 | 0.61 | 83.2 | 0.86 | 1.70 | 13.2 | 78.1 | 3.66 | 0 |
| 5 | g1r0_1_p2 | g0r0_299_p7 | 9.94 | 0.64 | 80.5 | 0.84 | 1.73 | 14.1 | 76.2 | 2.15 | 0 |
| 6 | g1r0_5_p2 | g1r0_399_p1 | 8.24 | 0.62 | 84.8 | 0.84 | 1.62 | 13.9 | 75.2 | 3.93 | 0 |
| 7 | g1r0_5_p5 | g1r0_421_p4 | 8.24 | 0.75 | 85.0 | 0.83 | 1.77 | 13.2 | 74.9 | 3.37 | 0 |
| 8 | g0r0_6_p5 | g1r0_399_p1 | 8.94 | 0.83 | 82.8 | 0.84 | 1.55 | 13.7 | 74.4 | 1.79 | 0 |
| 9 | g1r0_3_p4 | g0r0_199_p5 | 8.42 | 0.87 | 84.9 | 0.84 | 1.68 | 14.6 | 75.7 | 2.43 | 0 |
| 10 | g1r0_6_p3 | g1r0_421_p4 | 9.42 | 1.04 | 85.4 | 0.83 | 2.22 | 14.8 | 77.3 | 2.43 | 0 |
| 11 | g0r0_2_p7 | g0r0_299_p7 | 9.13 | 0.60 | 83.5 | 0.82 | 1.79 | 14.8 | 73.1 | 3.88 | 0 |
| 12 | g1r0_2_p6 | g0r0_199_p5 | 8.24 | 0.87 | 84.9 | 0.83 | 1.71 | 14.8 | 73.1 | 2.68 | 0 |
| 13 | g1r0_5_p2 | g0r0_199_p5 | 8.67 | 0.71 | 82.7 | 0.83 | 1.85 | 13.3 | 77.8 | 2.20 | 0 |
| 14 | g1r0_0_p5 | g1r0_421_p4 | 9.48 | 0.78 | 83.9 | 0.82 | 1.80 | 15.1 | 75.5 | 2.80 | 0 |
| 15 | g1r0_4_p3 | g1r0_421_p4 | 8.99 | 0.77 | 82.0 | 0.82 | 2.07 | 14.9 | 73.5 | 2.24 | 0 |
| 16 | g1r0_3_p3 | g1r0_298_p5 | 9.11 | 0.41 | 85.2 | 0.80 | 2.08 | 15.0 | 73.8 | 1.95 | 0 |
| 17 | g1r0_0_p4 | g1r0_421_p4 | 9.64 | 0.91 | 83.8 | 0.82 | 1.91 | 13.8 | 74.3 | 2.31 | 0 |
| 18 | g0r0_4_p4 | g1r0_399_p1 | 9.43 | 0.76 | 82.5 | 0.82 | 1.70 | 14.5 | 73.4 | 3.41 | 0 |
| 19 | g0r0_6_p0 | g0r0_299_p7 | 9.01 | 0.64 | 83.2 | 0.81 | 2.10 | 15.0 | 73.1 | 2.92 | 0 |
| 20 | g0r0_6_p0 | g1r0_399_p1 | 8.79 | 0.75 | 84.7 | 0.80 | 2.00 | 15.2 | 71.9 | 2.83 | 0 |

### 5.2 Epitope and hotspot engagement (AF3 complex, 4.5 Å heavy atom)

Hotspots in UniProt P43220 numbering. Every submitted design contacts all six. Distances are the closest binder heavy atom to that residue.

| # | design | parent lead | ECL2 /3 | ECL3 /3 | total /6 | min dist Å | epitope res | T298 | R299 | N300 | D372 | R380 | L384 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | g1r0_1_p4 | g0r0_199_p5 | 3 | 3 | 6 | 2.11 | 58 | 3.67 | 2.97 | 3.13 | 4.40 | 2.11 | 3.89 |
| 2 | g1r0_5_p4 | g0r0_299_p7 | 3 | 3 | 6 | 2.44 | 46 | 2.44 | 2.59 | 2.71 | 2.94 | 3.01 | 3.26 |
| 3 | g1r0_1_p3 | g1r0_399_p1 | 3 | 3 | 6 | 2.65 | 50 | 4.07 | 4.39 | 2.65 | 3.00 | 2.70 | 3.60 |
| 4 | g0r0_5_p7 | g0r0_299_p7 | 3 | 3 | 6 | 2.63 | 46 | 3.19 | 2.63 | 2.88 | 3.08 | 2.74 | 3.21 |
| 5 | g1r0_1_p2 | g0r0_299_p7 | 3 | 3 | 6 | 1.98 | 40 | 1.98 | 2.96 | 2.95 | 2.70 | 2.53 | 3.69 |
| 6 | g1r0_5_p2 | g1r0_399_p1 | 3 | 3 | 6 | 1.90 | 54 | 3.62 | 2.97 | 4.20 | 2.94 | 1.90 | 3.73 |
| 7 | g1r0_5_p5 | g1r0_421_p4 | 3 | 3 | 6 | 2.42 | 49 | 3.40 | 4.11 | 2.91 | 2.42 | 3.54 | 3.69 |
| 8 | g0r0_6_p5 | g1r0_399_p1 | 3 | 3 | 6 | 2.24 | 54 | 3.67 | 2.89 | 3.66 | 3.60 | 2.24 | 3.76 |
| 9 | g1r0_3_p4 | g0r0_199_p5 | 3 | 3 | 6 | 2.19 | 65 | 3.61 | 2.82 | 3.22 | 3.36 | 2.19 | 3.77 |
| 10 | g1r0_6_p3 | g1r0_421_p4 | 3 | 3 | 6 | 2.24 | 45 | 3.41 | 2.51 | 2.90 | 3.11 | 2.24 | 4.07 |
| 11 | g0r0_2_p7 | g0r0_299_p7 | 3 | 3 | 6 | 3.12 | 40 | 3.23 | 3.69 | 4.24 | 3.12 | 3.54 | 3.39 |
| 12 | g1r0_2_p6 | g0r0_199_p5 | 3 | 3 | 6 | 2.28 | 61 | 3.55 | 3.33 | 3.68 | 3.16 | 2.28 | 3.61 |
| 13 | g1r0_5_p2 | g0r0_199_p5 | 3 | 3 | 6 | 2.36 | 60 | 3.56 | 3.24 | 3.62 | 3.44 | 2.36 | 3.65 |
| 14 | g1r0_0_p5 | g1r0_421_p4 | 3 | 3 | 6 | 1.79 | 42 | 3.62 | 2.93 | 3.08 | 3.97 | 1.79 | 3.85 |
| 15 | g1r0_4_p3 | g1r0_421_p4 | 3 | 3 | 6 | 2.28 | 45 | 3.55 | 3.25 | 2.52 | 2.76 | 2.28 | 3.57 |
| 16 | g1r0_3_p3 | g1r0_298_p5 | 3 | 3 | 6 | 2.64 | 46 | 2.96 | 3.12 | 3.78 | 3.05 | 2.64 | 3.66 |
| 17 | g1r0_0_p4 | g1r0_421_p4 | 3 | 3 | 6 | 2.56 | 47 | 4.02 | 3.26 | 4.24 | 3.98 | 2.56 | 3.84 |
| 18 | g0r0_4_p4 | g1r0_399_p1 | 3 | 3 | 6 | 2.24 | 49 | 3.89 | 2.70 | 4.18 | 3.56 | 2.24 | 3.83 |
| 19 | g0r0_6_p0 | g0r0_299_p7 | 3 | 3 | 6 | 2.65 | 45 | 3.57 | 2.65 | 3.44 | 3.40 | 2.92 | 3.24 |
| 20 | g0r0_6_p0 | g1r0_399_p1 | 3 | 3 | 6 | 2.63 | 52 | 2.63 | 4.38 | 3.32 | 3.20 | 2.83 | 3.48 |

### 5.3 Interface energetics (Rosetta InterfaceAnalyzer)

Computed on the AF3 complex. ddG is also given per 1,000 Å² of buried surface, since a large raw ddG on an oversized interface reflects area rather than quality.

| # | design | parent lead | ddG | ddG / 1k Å² | dSASA_int Å² | shape compl. | packstat | H-bonds | interface res | hydrophobic % |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | g1r0_1_p4 | g0r0_199_p5 | -171.4 | -31.3 | 5481 | 0.62 | 0.58 | 27 | 48 | 56.2 |
| 2 | g1r0_5_p4 | g0r0_299_p7 | -111.4 | -27.7 | 4019 | 0.61 | 0.58 | 20 | 36 | 38.9 |
| 3 | g1r0_1_p3 | g1r0_399_p1 | -111.6 | -27.4 | 4065 | 0.65 | 0.56 | 13 | 38 | 47.4 |
| 4 | g0r0_5_p7 | g0r0_299_p7 | -94.8 | -27.2 | 3491 | 0.60 | 0.62 | 12 | 32 | 43.8 |
| 5 | g1r0_1_p2 | g0r0_299_p7 | -95.6 | -25.2 | 3798 | 0.64 | 0.65 | 9 | 35 | 45.7 |
| 6 | g1r0_5_p2 | g1r0_399_p1 | -117.8 | -28.1 | 4194 | 0.63 | 0.59 | 11 | 38 | 52.6 |
| 7 | g1r0_5_p5 | g1r0_421_p4 | -113.3 | -27.7 | 4095 | 0.64 | 0.53 | 9 | 39 | 64.1 |
| 8 | g0r0_6_p5 | g1r0_399_p1 | -132.1 | -24.6 | 5360 | 0.56 | 0.51 | 13 | 49 | 49.0 |
| 9 | g1r0_3_p4 | g0r0_199_p5 | -126.3 | -24.1 | 5237 | 0.59 | 0.56 | 21 | 49 | 51.0 |
| 10 | g1r0_6_p3 | g1r0_421_p4 | -134.3 | -26.2 | 5127 | 0.65 | 0.56 | 20 | 45 | 37.8 |
| 11 | g0r0_2_p7 | g0r0_299_p7 | -105.9 | -27.1 | 3915 | 0.64 | 0.55 | 15 | 34 | 44.1 |
| 12 | g1r0_2_p6 | g0r0_199_p5 | -122.2 | -26.4 | 4633 | 0.58 | 0.56 | 12 | 38 | 55.3 |
| 13 | g1r0_5_p2 | g0r0_199_p5 | -144.9 | -26.7 | 5421 | 0.60 | 0.54 | 12 | 46 | 50.0 |
| 14 | g1r0_0_p5 | g1r0_421_p4 | -114.8 | -27.5 | 4179 | 0.61 | 0.56 | 13 | 37 | 54.0 |
| 15 | g1r0_4_p3 | g1r0_421_p4 | -139.2 | -26.2 | 5309 | 0.58 | 0.56 | 19 | 49 | 51.0 |
| 16 | g1r0_3_p3 | g1r0_298_p5 | -118.2 | -29.5 | 4015 | 0.61 | 0.58 | 18 | 33 | 54.5 |
| 17 | g1r0_0_p4 | g1r0_421_p4 | -107.8 | -25.1 | 4291 | 0.61 | 0.56 | 9 | 40 | 57.5 |
| 18 | g0r0_4_p4 | g1r0_399_p1 | -111.0 | -24.3 | 4566 | 0.61 | 0.52 | 9 | 43 | 48.8 |
| 19 | g0r0_6_p0 | g0r0_299_p7 | -119.2 | -27.1 | 4391 | 0.60 | 0.55 | 18 | 43 | 44.2 |
| 20 | g0r0_6_p0 | g1r0_399_p1 | -105.2 | -25.8 | 4083 | 0.65 | 0.53 | 17 | 35 | 42.9 |

### 5.4 Designed sequence (ProteinMPNN) and developability

MPNN diversity is measured across the 8 sequences sampled for that design's backbone. No design carries a cysteine.

| # | design | parent lead | len | net charge pH 7.4 | pI | GRAVY | aromatic | Cys | N-glyc sequons | NG | MPNN distinct /8 | MPNN median id |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | g1r0_1_p4 | g0r0_199_p5 | 60 | -3.1 | 4.46 | 0.01 | 10.0% | 0 | 0 | 0 | 8 | 0.58 |
| 2 | g1r0_5_p4 | g0r0_299_p7 | 60 | -0.1 | 6.89 | -0.94 | 5.0% | 0 | 0 | 0 | 8 | 0.53 |
| 3 | g1r0_1_p3 | g1r0_399_p1 | 60 | -1.1 | 5.17 | -0.69 | 8.3% | 0 | 0 | 0 | 8 | 0.56 |
| 4 | g0r0_5_p7 | g0r0_299_p7 | 60 | -4.1 | 4.49 | -0.55 | 5.0% | 0 | 0 | 0 | 8 | 0.52 |
| 5 | g1r0_1_p2 | g0r0_299_p7 | 60 | -0.1 | 6.84 | -0.59 | 3.3% | 0 | 0 | 0 | 8 | 0.50 |
| 6 | g1r0_5_p2 | g1r0_399_p1 | 60 | -1.1 | 5.14 | -0.30 | 6.7% | 0 | 0 | 0 | 8 | 0.52 |
| 7 | g1r0_5_p5 | g1r0_421_p4 | 60 | +1.0 | 8.87 | -0.37 | 16.7% | 0 | 0 | 0 | 8 | 0.57 |
| 8 | g0r0_6_p5 | g1r0_399_p1 | 60 | +1.9 | 9.60 | -0.27 | 10.0% | 0 | 0 | 1 | 8 | 0.51 |
| 9 | g1r0_3_p4 | g0r0_199_p5 | 60 | -5.1 | 4.30 | -0.35 | 6.7% | 0 | 0 | 0 | 8 | 0.49 |
| 10 | g1r0_6_p3 | g1r0_421_p4 | 60 | -0.1 | 6.83 | -0.80 | 10.0% | 0 | 0 | 0 | 8 | 0.49 |
| 11 | g0r0_2_p7 | g0r0_299_p7 | 60 | -0.9 | 5.88 | -0.58 | 3.3% | 0 | 0 | 0 | 8 | 0.42 |
| 12 | g1r0_2_p6 | g0r0_199_p5 | 60 | -3.0 | 4.75 | -0.68 | 15.0% | 0 | 0 | 0 | 8 | 0.47 |
| 13 | g1r0_5_p2 | g0r0_199_p5 | 60 | +2.9 | 9.91 | -0.19 | 8.3% | 0 | 1 | 1 | 8 | 0.54 |
| 14 | g1r0_0_p5 | g1r0_421_p4 | 60 | -1.1 | 5.04 | -0.23 | 16.7% | 0 | 0 | 0 | 8 | 0.53 |
| 15 | g1r0_4_p3 | g1r0_421_p4 | 60 | -1.1 | 4.95 | -0.12 | 11.7% | 0 | 0 | 0 | 8 | 0.51 |
| 16 | g1r0_3_p3 | g1r0_298_p5 | 60 | -0.1 | 6.86 | -0.68 | 11.7% | 0 | 0 | 0 | 8 | 0.55 |
| 17 | g1r0_0_p4 | g1r0_421_p4 | 60 | +1.9 | 9.44 | -0.29 | 11.7% | 0 | 0 | 0 | 8 | 0.53 |
| 18 | g0r0_4_p4 | g1r0_399_p1 | 60 | -2.1 | 4.92 | -0.93 | 10.0% | 0 | 0 | 0 | 8 | 0.55 |
| 19 | g0r0_6_p0 | g0r0_299_p7 | 60 | -1.1 | 5.18 | -0.67 | 5.0% | 0 | 0 | 0 | 8 | 0.57 |
| 20 | g0r0_6_p0 | g1r0_399_p1 | 60 | +0.9 | 8.99 | -0.57 | 8.3% | 0 | 0 | 1 | 8 | 0.51 |

### 5.5 Sequences

| # | name | sequence |
|---|------|----------|
| 1 | glp1r_01_g0r0_199_p5_g1r0_1_p4 | `MGILEFTELLSLALEARGWRALAAQLAPAALADDESWAKYQAVLARSRAAERELERLFYE` |
| 2 | glp1r_02_g0r0_299_p7_g1r0_5_p4 | `GIGELTFTDLLASSLRAKAYKEMIKRLDQKAIENPEKEAEYRAKAEELRKEREKDLKVLE` |
| 3 | glp1r_03_g1r0_399_p1_g1r0_1_p3 | `SIGELTFTNLLALSYEVKALQARLERAKEELLRKGPEFEPEYLKVKKELEEARKKLDDAW` |
| 4 | glp1r_04_g0r0_299_p7_g0r0_5_p7 | `SIGELTFTDLYAAALEAKAYADQIKRLEEEAIRNPETAAELEAKAEELRKERARAKAILE` |
| 5 | glp1r_05_g0r0_299_p7_g1r0_1_p2 | `SIGVLTFTNLYSNVLRSKALRAQIRRLEQESLENPSLAPELEKRIEELRAERDRAESLLN` |
| 6 | glp1r_06_g1r0_399_p1_g1r0_5_p2 | `SLGELEFVDLYSLALEVKALKATVERLKKELLETGTEAEPKYIAAKKALEEAKKALESKW` |
| 7 | glp1r_07_g1r0_421_p4_g1r0_5_p5 | `GFGELSYTNLLYHALEAKAWTDFAKALLPKALADPAYEPQYWQAKEKAREAKAKLLSIWY` |
| 8 | glp1r_08_g1r0_399_p1_g0r0_6_p5 | `GLGELTFTDLYSLALEAKAAQRRVARYKELLLKNGDADFPKYLAAVAKLKALRQELEERF` |
| 9 | glp1r_09_g0r0_199_p5_g1r0_3_p4 | `MGELSFSELLSEALEARALKAQAEALKEAALESDEAAARYERAVAALREARDRLNELWYG` |
| 10 | glp1r_10_g1r0_421_p4_g1r0_6_p3 | `SIYESVFTELLGLSLEAKTLEKQLKTVEKLAEEDPSYQPKLDWLKKQLKEAKGRMESRWY` |
| 11 | glp1r_11_g0r0_299_p7_g0r0_2_p7 | `ALGELTFTDLLSNSLRAKALKEMIKRAEQAIIEHPELEPELLPRIEKARKEYKEAKSSVE` |
| 12 | glp1r_12_g0r0_199_p5_g1r0_2_p6 | `MGELTFTDLYADMLAALSYREFYEIQLKKYLEDPSNEHLREEAKARYQAARDALIKRFYE` |
| 13 | glp1r_13_g0r0_199_p5_g1r0_5_p2 | `MGLLEFTELLSLALEVKGLTKLLAIKAKEADKNGSDWPKYVAAKEQLKAAKKKLVDTWYG` |
| 14 | glp1r_14_g1r0_421_p4_g1r0_0_p5 | `SLSESFFTELYGLALEAASYAKFLAIIEKKYLEDPSFEPKYLEAKARLKAARDELLSRWY` |
| 15 | glp1r_15_g1r0_421_p4_g1r0_4_p3 | `AIGELTFTELLSLSLEVRSLTAFLKQFEPEYLANPSYEPKYLALKAQLDKAKTALLERWN` |
| 16 | glp1r_16_g1r0_298_p5_g1r0_3_p3 | `MGELEFTNYMAAVLGYKAYLTMAKRYREAAEFASPEKRAELLAKAEELEEKRKKYKEIIE` |
| 17 | glp1r_17_g1r0_421_p4_g1r0_0_p4 | `SLSEGFFTELLSLALEAKSLKETLAELKKKAEADPSYLPYYLKAKAALKEAKRKLLEAWY` |
| 18 | glp1r_18_g1r0_399_p1_g0r0_4_p4 | `SLGELTFTDLLALALEAKALKRRAERWREKYDLEGEEAEKKYLEAIEKYKEVKEKLDEKF` |
| 19 | glp1r_19_g0r0_299_p7_g0r0_6_p0 | `SIGELVFTEYLAAALEKKAAREAIKQLGEQAEKNPEDEAAVRAKMAEWRKRLKEAEKILE` |
| 20 | glp1r_20_g1r0_399_p1_g0r0_6_p0 | `GLGELTFTDLMYLALEKKALERRLARAEKLLLENGDKYEPLYLAVRKELRATREELERRF` |

<!-- END AUTOGENERATED LEAD TABLE -->
