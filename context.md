# Human GLP1R agonist design: challenge context and evidence-based strategy

Updated: 2026-10-17
Status: literature-grounded strategy; residue candidates are provisional until coordinate-level validation.

## 1. Purpose and scope

Develop de novo protein binders that both bind human GLP1R and activate its Gs/cAMP pathway. Affinity and agonism are separate objectives. The central design question is not simply where a protein can bind, but whether its binding geometry can productively engage receptor activation machinery. [C1, P1, P2]

This document consolidates the challenge specification and three uploaded primary papers. It replaces overconfident statements in the earlier discussion with evidence-ranked hypotheses. No receptor coordinates were supplied or analyzed in this session: no new distances, solvent-accessibility measurements, docking results, or state-selectivity energies are claimed. The existing SharePoint context.md and memory.md were not accessible and have not been merged. This is a new standalone context.md, not an edit of an inspected repository file.

Keep scientific context and design rationale here. Keep code versions, bugs, campaign parameters, candidate identifiers and experimental/computational outcomes in memory.md. Do not invent campaign history.

## 2. Challenge specification

The following is verified against the challenge page on 2026-10-17. [C1]

- Objective: design a human GLP1R agonist; binding without functional activation does not satisfy the challenge.
- Reference: UniProt P43220-1, residues 24–432, 409 amino acids.
- Purified target: supplier construct 13944-H18H-DD, N-terminal FLAG and C-terminal polyhistidine tags, in DDM/CHS detergent.
- Construct uncertainty: proprietary thermostabilizing mutations are not disclosed. The displayed sequence is the wild-type UniProt segment, not necessarily the exact assay-protein sequence.
- Organizer-recommended epitope: peptide-binding groove of the extracellular domain (ECD).
- Organizer structural reference: PDB 6X18, receptor author chain R.
- Binding is assessed against purified receptor in detergent. Confirmed binders proceed to a cell-based activation assay.
- Functional readout: human cells expressing surface GLP1R and a cAMP-response-element-controlled secreted alkaline phosphatase reporter; luminescence is compared with untreated cells and a reference agonist. This is an integrated downstream reporter, not a direct measurement of receptor conformation or Gs binding.
- Single-chain protein length: 10–250 amino acids. Accepted formats include proteins, nanobodies, scFvs and Fabs; format-specific rules apply.
- De novo designs only. Track 1: 20–40 submissions according to the FAQ; Tracks 2/3: at most 20.
- Deadline listed: 2026-10-25, 23:59 Anywhere on Earth. Recheck rules and deadline before submission.

Implication: a candidate must tolerate two different target environments. A detergent-compatible binder may fail on glycosylated cell-surface receptor; a membrane-dependent pose may fail the initial purified-target binding gate. [C1; design inference]

## 3. Evidence and terminology

Use these labels in future updates:

- **Published structural evidence:** contacts or conformations reported in a primary structure paper.
- **Published functional evidence:** measured expression, affinity, potency, efficacy or pathway bias.
- **Coordinate-derived evidence:** measurements performed on explicitly identified coordinate files; currently not available here.
- **Design hypothesis:** a proposed binding mode or mechanism requiring testing.

A native ligand-contact residue is not automatically an energetic hotspot. A loss-of-function receptor mutation is not automatically evidence that a binder should contact that residue. Mutation can alter folding, trafficking, receptor abundance, allosteric coupling or ligand recognition. [P2]

“Hotspot” below means a candidate receptor contact seed for design, not a validated high-success agonist epitope. No numerical success probability is established.

## 4. Structural target selection

| Structure | Role | Evidence and limitations |
|---|---|---|
| **6X18** | Primary peptide-like active-state template | Human GLP1R–GLP-1–Gs complex, 2.1 Å global resolution. Use receptor author chain R, verifying author/label chain mapping in mmCIF. The receptor sequence is not thermostabilized as in 6LN2, but the structural preparation includes engineered tags/signal sequence, dominant-negative Gs and Nb35. Global resolution does not imply uniform local accuracy. [P1, C1, S1] |
| **6LN2** | Inactive-state comparator | Peptide-free, closed-ECD human receptor, 3.2 Å, with Fab7F38, a negative allosteric modulator, rubredoxin insertion and 11 receptor mutations. Not an unmodified apo receptor or a definitive inactive ensemble. [P3] |
| **6X1A** | Alternative active mechanism | PF-06882961-bound receptor from the same study; useful for water-mediated activation and ECD/TM7 rearrangements. Not interchangeable with the peptide-bound geometry. [P1] |
| **6X19** | cAMP-active alternative geometry | CHU-128-bound receptor; useful because substantial cAMP activity occurs with different extracellular geometry and signaling bias. [P1] |
| **6B3J / 6ORV** | Optional active-state diversity | Exendin-P5 / TT-OAD2 comparators discussed in P1. Use only after inspecting their construct and local structural quality. [P1] |

Do not use the older rabbit structure 5VAI as the primary source of human atomic contact assignments. P1 explicitly revises several extracellular contacts relative to that lower-resolution model. [P1, Figure S4]

## 5. What activation requires

### 5.1 Canonical peptide mechanism

The GLP-1 C-terminal region anchors to the ECD; its N-terminal region engages the extracellular transmembrane cavity. Productive engagement supports receptor-core rearrangement, the TM6 kink and intracellular opening compatible with Gs coupling. ECD affinity alone does not demonstrate this coupling. [P1–P3]

For native GLP-1, the N-terminal H7/A8/E9 region occupies the deep pocket. H7 contacts the TM3/TM5 region; E9 contacts Y152 and R190 and participates in water-supported interactions. GLP-1 also engages ECL2 and ECL3/TM7. [P1, Figure 3 and Table S2]

### 5.2 Important exception

Direct physical mimicry of every GLP-1 N-terminal contact is **not universally required for agonism**. PF-06882961 activates through a shallower ligand pose plus structural waters and receptor rearrangements; CHU-128 also activates cAMP with a different pose and signaling profile. Thus deep-pocket penetration is a strong peptide-mimetic hypothesis, not an absolute requirement for every designed agonist. [P1]

Likewise, the fully extended ECD geometry in 6X18 is not the only agonist-compatible extracellular arrangement. Opening the ECD alone has not been shown here to be sufficient for GLP1R activation. [P1, P3]

### 5.3 Consequence for the challenge

Prioritize a credible Gs/cAMP activation mechanism. Full reproduction of GLP-1 arrestin recruitment, calcium signaling or internalization is not a stated requirement. However, cAMP reporter activity alone does not establish a complete pharmacological profile or therapeutic equivalence. [C1, P1]

## 6. Residue-level evidence map

All positions below use human receptor numbering as reported in the papers, not cropped design-file indices. GLP-1 positions use conventional hormone numbering, beginning at H7 for GLP-1(7–36). The contact assignments below are published observations from P1 Table S2, not new coordinate measurements. [P1]

### 6.1 ECD anchoring candidates

| Receptor residue | Published GLP-1 contact | Design interpretation |
|---|---|---|
| L32 | E21, A24, F28; hydrogen bond reported with E21 | Candidate at the N-terminal ECD helix; distinguish backbone and side-chain contributions before constraining. [P1, Table S2] |
| T35 | F28, I29 | Supporting anchoring contact. [P1, Table S2] |
| W39 | F28, L32, R36; aromatic interactions reported | Strong candidate to inspect for a native-groove anchoring campaign, conditional on accessible side-chain area. [P1, Table S2; hypothesis] |
| E68, Y69 | E68: L32/R36; Y69: I29/L32/V33 | Candidate supporting groove contacts; not proven energetic hotspots. [P1, Table S2] |
| P90, W91 | I29 | Candidate hydrophobic groove neighborhood. [P1, Table S2] |
| R121 | V33; hydrogen bond reported | Candidate polar anchoring contact; inspect the contacting atoms before designing chemistry. [P1, Table S2] |

**Provisional ECD seed pool:** W39, Y69, W91, R121. Choose small spatially compatible subsets only after coordinate inspection. Do not assume these four residues can all be contacted by a given scaffold. ECD-only binding is an affinity/control strategy, not a demonstrated agonist mechanism. [P1; hypothesis]

W33 is a different case: it is critical for the small molecules PF-06882961 and CHU-128, whereas W33A did not affect GLP-1-mediated cAMP in P1. Do not label W33 a universal peptide-agonism hotspot. [P1, Figure 7 and Figure S5]

### 6.2 ECL2 and adjacent TM5

| Residue | Structural/functional evidence | Design interpretation |
|---|---|---|
| **T298** | Contacts GLP-1 T13/S14/S17 and water-linked G10. T298A changes efficacy despite no significant agonist-affinity loss; increased apparent potency is confounded by increased functional receptor abundance. [P1 Table S2; P2] | Primary candidate at the peptide-facing loop surface, pending accessibility. |
| **R299** | Contacts S14/S17/S18/E21. R299A alters affinity, cAMP efficacy and pathway bias. [P1 Table S2; P2] | Primary polar-contact candidate; preserve a productive loop/network geometry rather than merely maximizing electrostatics. |
| **N300** | Contacts G10/S14. N300A reduces affinity and cAMP coupling despite increased antibody-detected surface expression. [P1 Table S2; P2] | Primary candidate with structural and functional support; not proof that external engagement causes activation. |
| W297 | P2 reports strong functional effects, structural conservation and prior cross-link proximity. It is not listed as a direct GLP-1 contact in P1 Table S2. [P1, P2] | Secondary/structural candidate, not a mandatory exposed hotspot. |
| N302, M303, N304, Y305, L307 | Functional effects on binding, coupling or bias; N304A particularly separates efficacy from affinity. [P2] | Mechanistic context and possible supporting residues; direct accessibility/contact suitability remains unverified. |
| **W306** | P1 places it in **TM5**, contacting H7/G10/T11. W306A eliminates detectable surface expression in P2. [P1, P2] | Deep pocket/structural residue, not a clean isolated activation switch or automatically accessible ECL2 hotspot. |
| C296 | Conserved disulfide linkage with C226; mutation affects structure/function. [P2] | Preserve receptor structural integrity. Do not prioritize the disulfide as a generic binder-contact seed. |

**Provisional ECL2 seed pool: T298, R299, N300.** This is a more defensible starting set than W297–R299 alone. Whether to use one, two or all three depends on the actual interface geometry. [P1, P2; hypothesis]

### 6.3 ECL3/TM7: activation-relevant contacts omitted from the earlier emphasis

| Residue | Evidence | Design interpretation |
|---|---|---|
| **R380** | GLP-1 D15 ionic interaction and T11 contact; PF-06882961 also uses an ionic interaction here. [P1 Table S2] | High-priority complementary activation-contact candidate; test an appropriately placed acidic group rather than assuming any salt bridge is agonistic. |
| **L384** | Contacts GLP-1 A8/T11/D15 and PF-06882961. [P1 Table S2] | Supporting hydrophobic candidate, conditional on access. |
| D372 | GLP-1 T11 contact; part of an activation-associated polar network. [P1] | Candidate supporting interaction/network constraint. |
| E373, K383, R310 | P1 discusses strong signaling effects and polar-network roles, not a universal requirement for direct ligand contact. [P1] | Monitor network compatibility; do not force all as binder hotspots. |
| E387, L388 | GLP-1 N-terminal pocket contacts. [P1 Table S2] | Deep-pocket geometry constraints for penetrating designs, not default exposed surface seeds. |

**Provisional complementary seed:** R380, with L384 and D372 as context-dependent support. An ECL2-plus-R380 design is a testable hypothesis, not a geometrically validated epitope. [P1; hypothesis]

### 6.4 Deep orthosteric pocket and upper TM contacts

- **Q234, V237, W306, I309, R310, I313:** H7-associated contact/network neighborhood. Table S2 adds I309 and identifies R310–H7 as water mediated. [P1]
- **Y152, R190:** direct E9 hydrogen-bond partners; **Y241** participates through waters. [P1]
- **L141, L144, Y148:** F12 hydrophobic/aromatic neighborhood. [P1]
- **K197:** T13 hydrogen bond; **Y205:** S17/L20/E21 contacts. [P1]
- **Q210, W214, Q221:** ECL1 contacts with more C-terminal GLP-1 positions; W214 contacts F28/W31 in Table S2. [P1]

These are mechanistic landmarks, not a single mandatory contact list. A globular binder may not reach deep-pocket residues without a narrow protrusion. [P1; geometric hypothesis]

## 7. Proposed design campaigns

### Campaign A: peptide-trajectory-inspired orthosteric agonist

**Preferred mechanistic starting hypothesis**, not a validated best-performing approach.

- Use 6X18 as the primary template.
- Explore compact scaffolds with an exposed, narrow activation element directed into the orthosteric vestibule.
- Candidate surface seeds: T298/R299/N300, with an independent branch exploring R380 and supporting L384/D372 contacts.
- Inspect whether a protruding element can reproduce selected N-terminal GLP-1 interactions without membrane or receptor clashes.
- Add ECD anchoring contacts only when the scaffold can span them with realistic geometry; do not force an ECD–ECL2 bridge before measuring distances.
- Main risks: pocket inaccessibility, nonproductive loop immobilization, antagonism, dependence on a free N-terminal chemical group that the scaffold does not reproduce, and incompatibility with assay tags or construct mutations.

Rationale derives from the native contact topology and activation-associated networks. No tested de novo scaffold is established by the supplied papers. [P1–P3; hypothesis]

### Campaign B: upper-vestibule/network agonist

- Explore ECL2 plus ECL3/TM7-facing engagement without insisting on deep H7-like penetration.
- Use 6X1A and 6X19 as mechanistic comparators to identify alternative activation-compatible networks, not as interchangeable fixed targets.
- Ask whether the designed interface supports productive polar packing and hydration rather than simply occupying a large surface.
- Main risk: small-molecule mechanisms may not transfer to a protein scaffold; docking cannot establish water-network stabilization or efficacy.

This branch is motivated by non-peptide agonism in P1. [P1; hypothesis]

### Campaign C: ECD anchoring/control branch

- Inspect W39/Y69/W91/R121 as provisional native-groove seed candidates.
- Useful for affinity-focused designs and for testing whether ECD engagement alone produces activity.
- Do not describe an ECD-only binder as an agonist before functional data.
- A conformational clamp spanning ECD and extracellular loops remains exploratory. P3 shows that locking a closed arrangement impairs GLP-1 potency; it does not prove that locking an open arrangement is sufficient for activation.

[P1, P3; hypothesis]

No evidence-based numerical allocation or success rate can yet be assigned across campaigns.

## 8. State comparison: useful, but not an agonism classifier

Compare poses with 6LN2 after explicit sequence/numbering reconciliation and accounting for engineered mutations and fusion insertion. Align the transmembrane core separately from the ECD; whole-receptor alignment can obscure domain motion. [P3; proposed analysis]

Do not automatically reject a candidate because it clashes with rigid 6LN2: an agonist can bind through conformational selection or induced fit. Conversely, a favorable active-state docking score does not demonstrate that the ligand shifts the receptor equilibrium or activates Gs. Gs and structural stabilizers already help maintain the active structures. [P1, P3]

Any active-minus-inactive score must specify score direction, target preparation and normalization. It is a comparative modeling heuristic, not a measured free-energy difference or calibrated probability of agonism.

## 9. Required coordinate-level work before production hotspot constraints

1. Obtain full mmCIF files for 6X18 and 6LN2; preferably add 6X1A and 6X19.
2. Map UniProt positions to author residue IDs, label residue IDs and cropped design indices. Explicitly handle the rubredoxin insertion and renumbering in 6LN2.
3. Check residue identities, missing atoms/residues, alternate conformations, local quality and construct modifications.
4. Recompute GLP-1 contacts and compare with P1 Table S2; distinguish direct heavy-atom contacts, hydrogen-bond geometry and water-mediated networks.
5. Measure solvent-accessible side-chain area and spatial distances among candidate seed sets.
6. Check membrane orientation, detergent-exposed hydrophobic surfaces, glycans and tag proximity. Absence of a modeled glycan does not prove absence on the cell-surface receptor.
7. Evaluate candidate approach vectors and clashes using the full receptor, not just a cropped epitope.
8. Export an auditable hotspot table with receptor position, residue identity, coordinate ID, contact atoms, accessibility, evidence source, uncertainty and campaign role.

Suggested contact cutoffs should be reported as analysis choices and sensitivity-tested; do not silently equate a distance threshold with an energetic hotspot. P1 Table S2 reports receptor residues within 4 Å of ligand, with additional interaction annotations. [P1]

## 10. Candidate evaluation and functional interpretation

Use independent axes rather than a single docking rank:

- Foldability and expression/developability.
- Interface plausibility and accessible epitope engagement.
- Activation-mechanism plausibility.
- Robustness across receptor conformations and assay environments.
- Diversity of scaffold, pose and mechanism.

Structural prediction confidence is not affinity, and neither confidence nor affinity is efficacy. Prioritize experimental functional testing rather than claiming agonism computationally. [P1, P2; design principle]

For interpretation of functional hits, distinguish receptor-dependent activation from reporter interference or nonspecific cellular effects. Compare concentration-dependent activity and maximal response with the reference agonist; use orthogonal cAMP evidence and receptor-negative controls where feasible. A ligand that only enhances a reference agonist response is not necessarily an autonomous agonist. These are proposed validation principles, not a claim about additional organizer assays. [C1; proposed validation]

## 11. Corrections to the earlier discussion

1. Do not transfer the older 5VAI contact model directly to 6X18. In the higher-resolution human structure, H7 occupies the TM3/TM5 neighborhood, and E9 directly engages Y152/R190; the earlier H7–R299 and Y145 emphasis is not the appropriate primary 6X18 assignment. [P1, Figure S4 and Table S2]
2. W306A loss of signaling is confounded by absent surface expression. It is not clean evidence for a directly targetable activation switch. [P2]
3. W297 is functionally important but is not a listed direct GLP-1 contact in P1 Table S2. [P1, P2]
4. ECL3/TM7, especially R380 and its associated network, deserves explicit consideration alongside ECL2. [P1]
5. An ECD–ECL2 bridge is a proposed architecture, not a validated high-success epitope.
6. Deep peptide-like penetration and a fully extended ECD are not universal requirements for all GLP1R agonists. [P1]
7. Use 6X1A and 6X19 as the directly relevant alternative agonist structures from the supplied study; 7LCI is not needed as the first comparator. [P1]
8. No quantitative success-rate ranking is justified by these sources.

## 12. Current decision summary

- **Primary structure:** 6X18, receptor author chain R; verify coordinate mapping before use. [C1, P1]
- **Inactive comparator:** 6LN2, with engineering and conformational-ensemble caveats. [P3]
- **First surface seed pool to validate:** T298, R299, N300. [P1, P2; hypothesis]
- **Complementary activation-facing candidate:** R380; inspect L384 and D372 as support. [P1; hypothesis]
- **ECD anchoring pool:** W39, Y69, W91, R121, conditional on geometry and accessibility. [P1; hypothesis]
- **Deep-pocket landmarks:** Q234/W306 neighborhood and Y152/R190/E387/L388 region, for penetrating designs rather than default exposed seeds. [P1]
- **Not yet justified:** compulsory hotspot sets, exact tool-specific indices, interface distances, high-success claims, or an assertion that any proposed binder is agonistic.

## 13. Sources and provenance

**[C1] Challenge specification.** Proteinbase, Anthropic × Adaptyv 2026, Challenge 03: GLP1R. Accessed 2026-10-17.
https://proteinbase.com/competitions/anthropic-adaptyv-2026/challenges/glp1r

**[S1] Structure record.** RCSB PDB 6X18. Accessed 2026-10-17.
https://www.rcsb.org/structure/6X18
Structure DOI: 10.2210/pdb6X18/pdb.

**[P1] Zhang et al.** Differential GLP-1R Binding and Activation by Peptide and Non-peptide Agonists. Molecular Cell 80, 485–500.e7 (2020).
[DOI:10.1016/j.molcel.2020.09.020]
User-uploaded file: mmc6.pdf. Despite its filename, this contains the main article, methods and supplementary material. Particularly relevant: main article pp. 489–495; Table S2, PDF pp. 28–30; Figures S4/S5. Contact assignments in this document are literature-derived from this source.

**[P2] Koole et al.** Second Extracellular Loop of Human Glucagon-like Peptide-1 Receptor (GLP-1R) Has a Critical Role in GLP-1 Peptide Binding and Receptor Activation. Journal of Biological Chemistry 287, 3642–3658 (2012).
[DOI:10.1074/jbc.M111.309328]
User-uploaded file: zbc3642.pdf. Particularly relevant: surface-expression results; Tables 1/2; discussion of W306A, C296 and indirect agonist-affinity effects. Some numerical table content is image-based in the extraction; exact numerical effect sizes have not been independently transcribed here.

**[P3] Wu et al.** Full-length human GLP-1 receptor structure without orthosteric ligands. Nature Communications 11, 1272 (2020).
[DOI:10.1038/s41467-020-14934-5]
User-uploaded file: s41467-020-14934-5.pdf. Particularly relevant: Figures 2/3, inactive construct engineering, disulfide-lock experiments and discussion of ECD flexibility.

## 14. Update policy

Preserve provenance. Add coordinate-derived measurements only with file identifiers, mapping conventions and analysis settings. Replace provisional seeds only when evidence warrants it. Record campaign-specific successes, failures and code changes in memory.md and link the relevant finding back here. Never convert a computationally plausible mechanism into a claim of experimentally demonstrated agonism.
