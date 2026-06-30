# CRISPR-Cas9 Knockout Design for ABL1 — with genome-wide off-target validation

A complete computational design of a CRISPR-Cas9 knockout strategy for the human **ABL1** gene (the leukemia driver, c-Abl), built from the gene sequence up: guide RNA discovery, on-target efficiency scoring, CFD off-target analysis, **genome-wide validation against CRISPOR**, and a wet-lab validation protocol with PCR primers.

This is the third in a set of projects centered on the same target. The companion projects [dock drugs into c-Abl](https://github.com/pitamberdasla108-ui/gpu-virtual-screening) and [predict c-Abl binding with deep learning](https://github.com/pitamberdasla108-ui/kinase-binding-predictor); here the same kinase is approached from the genome side — designing how to delete it.

---

## What this project does

Designs single-guide RNAs (sgRNAs) to knock out ABL1 by targeting an early coding region, then **validates the design the way a real lab would**: scoring guides for cutting efficiency, searching the whole genome for off-target sites, and designing the PCR assay to confirm the edit at the bench.

The headline result is a concrete demonstration of *why genome-wide off-target validation is non-negotiable*: the guides that scored best on sequence features alone turned out to be the **most dangerous** genome-wide — including one predicted to also cut an exon of *ABL2*, ABL1's close paralog.

---

## Approach & pipeline

**1. Fetch the real gene.** ABL1 reference transcript NM_005157.6 pulled from NCBI (5,578 bp mRNA; CDS 3,393 bp → 1,130 aa — matching the c-Abl protein used in the companion projects).

**2. Choose the knockout target.** Exons 2-3 (early, fully-coding, before the kinase domain at ~residue 242). A frameshift here disrupts the protein before any catalytic function — a clean knockout.

**3. Find candidate guides** (`find_guides.py`). Scan both strands for the SpCas9 motif (20 nt protospacer + NGG PAM). 67 raw sites → 53 after quality filters (GC 40-75%, no poly-T Pol III terminator).

**4. Score on-target efficiency** (`score_guides.py`). Interpretable, literature-based features (GC optimum, seed-region composition, 5′-G, homopolymer avoidance, PAM-proximal purine). *Note: an initial version saturated (all top guides = 1.0); it was redesigned into a discriminating points-based scorer spreading scores 0.54-1.0.*

**5. CFD off-target scoring** (`cfd_offtarget.py`). Implements the Cutting Frequency Determination algorithm (Doench 2016): off-target activity = product of position-dependent mismatch-tolerance weights. Demonstrates the core CRISPR-specificity principle — PAM-distal mismatches are tolerated (guide still cuts), seed-region mismatches abolish cutting.

**6. Genome-wide validation against CRISPOR** (`crispor_validation.py`). The top guides were submitted to CRISPOR (gold-standard tool, exact CFD matrix + real hg19 genome search). **This is the critical step** — see findings below.

**7. Validation assay design** (`design_primers.py`). PCR primers flanking the cut site for the selected guide, designed for a T7E1 mismatch-cleavage assay.

---

## Key finding: sequence-based scoring is not enough

| Guide | On-target | CRISPOR MIT spec. | Off-targets | Verdict |
|---|---|---|---|---|
| GTCTGAGTGAAGCCGCTCGT | 0.88 | **98** | 27 | ✅ SAFE |
| AGCCGCTCGTTGGAACTCCA | 0.94 | **90** | 41 | ✅ SAFE |
| CATCACGCCAGTCAACAGTC | 0.89 | **90** | 55 | ✅ SAFE |
| GCTGCTCAGCAGATACTCAG | 1.00 | 65 | 227 | ⚠️ MARGINAL |
| GCTGAGTATCTGCTGAGCAG | 1.00 | 53 | 244 | ❌ RISKY |
| GAGAAACACTCCTGGTACCA | 1.00 | 42 | 164 | ❌ RISKY — hits **ABL2 exon** |

The three guides my on-target scorer ranked highest (1.00) were the **worst** on genome-wide specificity. The top one is predicted to also cut an **exon of ABL2** — ABL1's closest paralog, exactly the gene you must not co-disrupt in a clean ABL1 knockout. Sequence-based scoring is blind to this because it never searches the genome; only the genome-wide search (CRISPOR) reveals it. **Final guide selection therefore required both: on-target efficiency AND validated genome-wide specificity.**

---

## Final design (orderable)

- **Target:** ABL1 exons 2-3 (early coding, upstream of the kinase domain)
- **Selected guide:** `GTCTGAGTGAAGCCGCTCGT` + NGG PAM (CRISPOR MIT specificity 98; 27 off-targets, none in coding exons of concern)
- **Validation primers:**
  - Forward `5'-TGTTGGAGATCTGCCTGAAGCTGGT-3'` (Tm 59.9 °C)
  - Reverse `5'-ATGTAGTTGCTTGGGACCCAGCCTT-3'` (Tm 60.1 °C)
- **Assay:** 346 bp amplicon; T7E1 mismatch cleavage yields ~137 + ~209 bp fragments in edited cells (uncut control stays a single 346 bp band).

---

## Repository contents
## Honest scope & limitations

- Guide design uses the canonical mRNA transcript; intron/exon boundaries at the genomic level would be cross-checked for a real experiment.
- The CFD implementation uses the real algorithm with a position-weight model; CRISPOR applies the exact published 400-value matrix + genome search (which is why it's the validation authority here).
- Validation primers are designed on transcript coordinates for demonstration; a real assay places them on genomic DNA (introns included) and BLAT-checks them.
- This is a **design and computational-validation** project — the wet-lab execution (transfection, T7E1, sequencing) is specified but not performed.

## References
- Doench et al., *Optimized sgRNA design to maximize activity and minimize off-target effects*, Nat Biotechnol 2016 (CFD).
- Haeussler et al., *Evaluation of off-target and on-target scoring algorithms and integration into CRISPOR*, Genome Biol 2016.
- Hsu et al., *DNA targeting specificity of RNA-guided Cas9 nucleases*, Nat Biotechnol 2013 (MIT score).
