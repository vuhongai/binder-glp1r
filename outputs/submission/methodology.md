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