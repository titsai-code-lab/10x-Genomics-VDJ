# 10x Genomics V(D)J IgG Cloning — Multi-Species Interactive Dashboard

An interactive **Plotly Dash** reference dashboard for single-cell B cell receptor (BCR) V(D)J sequencing using the **10x Genomics Chromium Next GEM Single Cell 5' Reagent Kits v2**. Covers primer sets, sequencing procedure, gene architecture maps, and isotype coverage for **mouse**, **human**, and **rabbit** IgG cloning.

> **Document reference:** CG000330 Rev A (August 2020) — *Chromium Next GEM Single Cell 5' Reagent Kits v2 (Dual Index)*
> **For Research Use Only. Not for diagnostic procedures.**

---

## Quick Start

### Requirements

```
Python ≥ 3.8
dash
plotly
pandas
```

### Installation & Launch

```bash
pip install dash plotly pandas
python vdj_mouse_bcr_dashboard.py
```

Then open **http://127.0.0.1:8050** in your browser.

---

## Dashboard Overview

The app has six tabs, each covering a distinct aspect of the workflow.

### 📋 Procedure
Step-by-step protocol cards for all **8 workflow stages**, from single-cell loading through Illumina sequencing and Cell Ranger analysis. Each card shows:
- Input materials and reagents
- PCR conditions and cycle numbers
- Cleanup and QC checkpoints
- Expected output size/peak

An interactive timeline chart provides an at-a-glance overview of the full workflow.

### 🧪 Primer Sequences
A filterable, sortable table of all primers with a **species toggle** (see below). Displays:

| Column | Description |
|--------|-------------|
| Name | Primer identifier |
| Part Number | 10x Genomics PN or custom designation |
| Species | Universal / Mouse / Human / Rabbit |
| Specificity | Universal / Mouse-Specific / Human-Specific / Rabbit-Specific |
| Isotype | IgG (heavy) / IgK (light) / IgL (light) / — |
| Stage | Where in the workflow this primer is used |
| Role | Functional description |
| Sequence (5'→3') | Oligonucleotide sequence |
| Length | Primer length in nucleotides |

**Row color coding:**
- Normal text → Universal primers
- **Bold orange underline** → Mouse-specific primers
- **Bold blue underline** → Human-specific primers
- ***Bold purple italic underline*** → Rabbit-specific primers *(computationally designed; require experimental validation)*

**Isotype column color coding:**
- Green → IgG heavy chain primer
- Blue → IgK light chain primer
- Purple → IgL light chain primer

### 🎯 Isotype Coverage
Sunburst chart of mouse BCR isotype hierarchy (IGH / IGK / IGL) with a checklist table confirming which isotypes are covered by Mix 1 v2 and Mix 2 v2.

### 📦 Kit Components
Complete parts list of required 10x Genomics kits, reagents, and accessories with storage temperatures color-coded (−80°C / −20°C / 4°C / Ambient).

### 📊 Analytics
- Bar chart: number of primers per workflow stage
- Sankey diagram: nested PCR BCR signal enrichment vs background removal
- Sequencing specifications table comparing V(D)J and 5' GEX library parameters

### 🧬 Gene Map
Cartoon domain architecture viewer showing where Mix 1 (outer) and Mix 2 (inner) primers anneal relative to the immunoglobulin constant region domains. Select an isotype from the dropdown to view its domain layout, primer positions, and a highlighted representative sequence panel.

---

## Species Toggle (Primer Sequences Tab)

When the **🧪 Primer Sequences** tab is active, a toggle strip appears below the tab bar:

| Option | Shows |
|--------|-------|
| 🐭 Mouse + Universal | 9 universal + 6 mouse IgG primers |
| 🧑 Human + Universal | 9 universal + 6 human IgG primers |
| 🐇 Rabbit + Universal | 9 universal + 6 rabbit IgG primers |
| ⚗️ Universal only | 9 universal primers only |
| 📋 All primers | All 27 primers across all species |

---

## Primer Sets — IgG Cloning Only

This dashboard is configured for **IgG heavy chain + IgK/IgL light chain cloning only**. IgM, IgD, IgA, and IgE primers have been removed. Each species has 6 species-specific primers (3 outer Mix 1, 3 inner Mix 2):

| Chain | Mouse (PN-1000255) | Human (PN-1000253) | Rabbit (custom) |
|-------|-------------------|-------------------|-----------------|
| IgG heavy (outer) | Mix 1 v2 Cγ1/2a/2b/2c/3 | Mix 1 v2 Cγ1/2/3/4 | Rab-Mix1 Cγ CH2 |
| IgG heavy (inner) | Mix 2 v2 Cγ CH1 | Mix 2 v2 Cγ CH1 | Rab-Mix2 Cγ CH1 |
| IgK light (outer) | Mix 1 v2 Cκ | Mix 1 v2 Cκ | Rab-Mix1 Cκ1/2 |
| IgK light (inner) | Mix 2 v2 Cκ | Mix 2 v2 Cκ | Rab-Mix2 Cκ1/2 |
| IgL light (outer) | Mix 1 v2 Cλ | Mix 1 v2 Cλ | Rab-Mix1 Cλ |
| IgL light (inner) | Mix 2 v2 Cλ | Mix 2 v2 Cλ | Rab-Mix2 Cλ |

> IgK and IgL light chain primers are **retained** because every IgG antibody is paired with either a κ or λ light chain. Both are needed to recover complete paired heavy/light chain sequences from single B cells.

---

## Nested PCR Strategy

The V(D)J enrichment uses **two rounds of PCR**:

```
Full-length cDNA
    ↓  Round 1 PCR — Mix 1 v2 (outer reverse primers)
       Outer constant-region primers + universal forward primer
       → Enriched BCR amplicons (~600–1200 bp)
    ↓  Round 2 PCR — Mix 2 v2 (inner reverse primers, nested)
       Inner primers positioned 5' of Mix 1 within constant region
       → Nested-enriched amplicons (~600 bp) + P5 added
    ↓  Library construction → Illumina sequencing
```

The **universal forward primer** (`5'-GATCTACACTCTTTCCCTACACGACGC-3'`) is identical in both rounds and for all three species.

---

## Sequencing Specifications

| Parameter | V(D)J Library | 5' GEX Library |
|-----------|--------------|----------------|
| Library type | Paired-end | Paired-end |
| Recommended read length | 150 bp × 2 | 26×10×10×90 bp |
| Minimum reads per cell | ≥ 5,000 | ≥ 20,000 |
| Read 1 encodes | 16-bp barcode + 10-bp UMI + V(D)J 5' end | 16-bp barcode + 10-bp UMI |
| Read 2 encodes | Internal V(D)J fragment | cDNA insert |
| Index read | i7 (+ i5 dual) | i7 (+ i5 dual) |
| Analysis tool | `cellranger vdj` | `cellranger count` |

---

## Rabbit Primer Design Notes

Rabbit primers are **not commercially available** from 10x Genomics. The primers in this dashboard are computationally designed and require experimental validation before use.

### Design Basis

| Isotype | Source |
|---------|--------|
| IgG outer (Cγ CH2) | Adapted from RIGHC2, Lavinder et al. 2014 *PLoS ONE* 9:e101322; GenBank L29172 |
| IgG inner (Cγ CH1) | Adapted from RIGHC1, Lavinder et al. 2014 *PLoS ONE*; GenBank L29172 |
| IgK outer/inner (Cκ1/2) | Adapted from RIGκC design, Lavinder et al. 2014; GenBank AY550529 |
| IgL outer/inner (Cλ) | Designed from rabbit IGLC sequence; GenBank X04050 |

### Key Rabbit Immunoglobulin Biology

- **No IgD isotype** (unlike mouse and human)
- **Single IgG isotype** (no subclasses)
- **13 IgA isotypes** (not included in this IgG-only workflow)
- **IgK dominant**: κ light chains comprise >90% of rabbit light chains
- All rabbit primers are marked **‡** and displayed in purple italic in the table

### Recommended Validation Steps

1. BLAST each primer sequence against *Oryctolagus cuniculus* genome assembly **OryCun2.0** (NCBI) to confirm unique mapping to the target constant region
2. Test in a conventional RT-PCR assay using rabbit PBMC RNA before 10x library preparation
3. Verify amplicon size on Bioanalyzer (expected ~600 bp nested product)
4. Confirm BCR enrichment efficiency with a pilot 10x run at low cell input

---

## Workflow — 8 Steps

| Step | Name | Key Action |
|------|------|------------|
| 1 | GEM Generation & Barcoding | Single cells partitioned into GEMs; barcoded TSO captures 5' end of mRNA |
| 2 | Post GEM-RT Cleanup & cDNA Amplification | GEMs broken; Dynabeads cleanup; cDNA PCR amplification; SPRIselect QC |
| 3 | V(D)J PCR Round 1 | Mix 1 v2 outer reverse primers + universal forward; SPRIselect cleanup |
| 4 | V(D)J PCR Round 2 | Mix 2 v2 inner (nested) reverse primers; adds P5; SPRIselect + Bioanalyzer QC |
| 5 | V(D)J Library Construction | Fragmentation, end repair, TruSeq adapter ligation, dual index PCR |
| 6 | 5' GEX Library Construction | Parallel processing of remaining cDNA into gene expression library |
| 7 | Illumina Sequencing | Paired-end 150 bp × 2; ≥5,000 reads/cell recommended for V(D)J |
| 8 | Cell Ranger vdj Analysis | Barcode/UMI filtering, V(D)J assembly, CDR3 annotation, clonotype calling |

---

## Key Kit Part Numbers

| Kit | Part Number | Storage |
|-----|-------------|---------|
| Chromium Next GEM Single Cell 5' Kit v2 (16 rxns) | PN-1000263 | −20°C |
| Chromium Next GEM Single Cell 5' Gel Bead Kit v2 | PN-1000264 | −80°C |
| Library Construction Kit | PN-1000190 | −20°C |
| **Mouse BCR Amplification Kit** | **PN-1000255** | −20°C |
| Mouse B Cell Mix 1 v2 | PN-2000258 | −20°C |
| Mouse B Cell Mix 2 v2 | PN-2000259 | −20°C |
| **Human BCR Amplification Kit** | **PN-1000253** | −20°C |
| Human B Cell Mix 1 v2 | PN-2000254 | −20°C |
| Human B Cell Mix 2 v2 | PN-2000255 | −20°C |
| Chromium Next GEM Chip K Single Cell Kit (48 rxns) | PN-1000286 | Ambient |
| Dual Index Kit TT Set A (96 rxns) | PN-1000215 | −20°C |
| Dual Index Kit TN Set A (96 rxns) | PN-1000250 | −20°C |
| Dynabeads MyOne SILANE | PN-2000048 | 4°C |
| Poly-dT RT Primer | PN-2000007 | −20°C |

---

## References

1. **10x Genomics User Guide** — *Chromium Next GEM Single Cell 5' Reagent Kits v2 (Dual Index) with Feature Barcode technology for Cell Surface Protein & Immune Receptor Mapping*. Document CG000330 Rev A, August 2020.

2. **Lavinder JJ, Hoi KH, Reddy ST, Wine Y, Georgiou G** (2014). Systematic Characterization and Comparative Analysis of the Rabbit Immunoglobulin Repertoire. *PLoS ONE* 9(6):e101322. doi:10.1371/journal.pone.0101322
   — Source of rabbit IgG primers RIGHC1 (inner) and RIGHC2 (outer), and RIGκC light chain primer design.

3. **Esteves PJ et al.** (2019). Genetic Diversity of IGHM and IGHE in the Leporids. *Animals* 9(11):955. doi:10.3390/ani9110955
   — Source of rabbit IgE primers FE12/RE21 used as design templates.

4. **Gertz EM et al.** (2013). Accuracy and coverage assessment of *Oryctolagus cuniculus* genes encoding immunoglobulins in the whole genome sequence assembly (OryCun2.0). *Immunogenetics* 65:749–762.
   — Reference genome for primer BLAST validation.

5. **GenBank accession numbers used for rabbit primer design:**
   - L29172 — *O. cuniculus* IGHG constant region
   - J00666 — *O. cuniculus* IGHM constant region
   - AY386696 / AY386696.1 — *O. cuniculus* IGH locus (IgG, IgM, IgE, IgA)
   - AY550529 — *O. cuniculus* IGKC (kappa light chain constant region)
   - X04050 — *O. cuniculus* IGLC (lambda light chain constant region)

---

## File Structure

```
vdj_mouse_bcr_dashboard.py   # Main application (single file)
README.md                    # This file
```

---

## Disclaimer

Rabbit primer sequences are **computationally designed** based on published constant-region sequences and design principles from the literature. They have **not been experimentally validated** and are provided for research planning purposes only. All primers should be validated in your own laboratory before use in 10x Genomics library preparation.

Mouse and human primer sequences follow the design logic of the 10x Genomics CG000330 user guide. The exact proprietary sequences for PN-2000254, PN-2000255, PN-2000258, and PN-2000259 are listed in the Oligonucleotide Sequences Appendix (p. 79) of CG000330. The sequences shown in this dashboard are representative and based on published literature; they may differ from the exact 10x Genomics commercial sequences.
# 10x-Genomics-VDJ
