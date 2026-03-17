"""
10x Genomics Chromium Next GEM Single Cell V(D)J Reagent Kits v2
Mouse B Cell (BCR) — Interactive Dashboard
CG000330 Rev A (August 2020)

Run:
    pip install dash plotly pandas
    python vdj_mouse_bcr_dashboard.py
Then open http://127.0.0.1:8050 in your browser.
"""

import dash
from dash import dcc, html, dash_table
import plotly.graph_objects as go
import plotly.express as px
import pandas as pd

# ──────────────────────────────────────────────────────────────────────────────
# DATA
# ──────────────────────────────────────────────────────────────────────────────

PRIMERS = pd.DataFrame([
    # ── Universal / shared primers ─────────────────────────────────────────
    {
        "Name": "Gel Bead TSO (v1 & v2)",
        "Part Number": "PN-1000264 / PN-1000267",
        "Species": "Universal",
        "Specificity": "Universal",
        "Role": "Captures 5' end; adds cell barcode + UMI",
        "Sequence (5'→3')": "CTACACGACGCTCTTCCGATCT[16-bp BC][10-bp UMI]TTTCTTATATrGrGrG",
        "Length": "~55 nt + barcode/UMI",
        "Stage": "GEM-RT",
    },
    {
        "Name": "Poly-dT RT Primer",
        "Part Number": "PN-2000007",
        "Species": "Universal",
        "Specificity": "Universal",
        "Role": "Reverse transcription of poly-A mRNA",
        "Sequence (5'→3')": "AAGCAGTGGTATCAACGCAGAGTACTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTVN",
        "Length": "58 nt",
        "Stage": "GEM-RT",
    },
    {
        "Name": "Universal Forward Primer (Mix 1 & 2)",
        "Part Number": "—",
        "Species": "Universal",
        "Specificity": "Universal",
        "Role": "Shared forward primer for all V(D)J PCR rounds (both species)",
        "Sequence (5'→3')": "GATCTACACTCTTTCCCTACACGACGC",
        "Length": "27 nt",
        "Stage": "V(D)J PCR 1 & 2",
    },
    {
        "Name": "TruSeq Adapter (sense strand)",
        "Part Number": "PN-220026",
        "Species": "Universal",
        "Specificity": "Universal",
        "Role": "Ligation adapter — enables Illumina library construction",
        "Sequence (5'→3')": "GATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
        "Length": "33 nt",
        "Stage": "Library construction",
    },
    {
        "Name": "Illumina TruSeq Read 1 Primer",
        "Part Number": "—",
        "Species": "Universal",
        "Specificity": "Universal",
        "Role": "Sequencing primer — reads barcode, UMI, and V(D)J 5' end",
        "Sequence (5'→3')": "ACACTCTTTCCCTACACGACGCTCTTCCGATCT",
        "Length": "33 nt",
        "Stage": "Sequencing",
    },
    {
        "Name": "Illumina TruSeq Read 2 Primer",
        "Part Number": "—",
        "Species": "Universal",
        "Specificity": "Universal",
        "Role": "Sequencing primer — reads internal V(D)J fragment",
        "Sequence (5'→3')": "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT",
        "Length": "34 nt",
        "Stage": "Sequencing",
    },
    {
        "Name": "TruSeq i5 Index Sequencing Primer",
        "Part Number": "—",
        "Species": "Universal",
        "Specificity": "Universal",
        "Role": "Index 2 read — sample demultiplexing",
        "Sequence (5'→3')": "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
        "Length": "33 nt",
        "Stage": "Sequencing",
    },
    {
        "Name": "Dual Index Kit TT — Forward Primer",
        "Part Number": "PN-3000431",
        "Species": "Universal",
        "Specificity": "Universal",
        "Role": "Adds P5 + i5 for Illumina cluster generation",
        "Sequence (5'→3')": "AATGATACGGCGACCACCGAGATCTACAC[10-bp i5]ACACTCTTTCCCTACACGACGCTC",
        "Length": "~63 nt",
        "Stage": "Library construction",
    },
    {
        "Name": "Dual Index Kit TT — Reverse Primer",
        "Part Number": "PN-3000431",
        "Species": "Universal",
        "Specificity": "Universal",
        "Role": "Adds P7 + i7 for Illumina cluster generation",
        "Sequence (5'→3')": "CAAGCAGAAGACGGCATACGAGAT[10-bp i7]GTGACTGGAGTTCAGACGTGT",
        "Length": "~55 nt",
        "Stage": "Library construction",
    },
    # ── Mouse B Cell specific primers (PN-1000255) ─────────────────────────
    {
        "Name": "Mouse B Cell Mix 1 v2 — IgM outer (Cμ)",
        "Part Number": "PN-2000258",
        "Species": "Mouse",
        "Specificity": "Mouse-Specific",
        "Role": "Round 1 outer reverse primer — anneals ~3' of CH1 in mouse Cμ constant region",
        "Sequence (5'→3')": "CAGGGGGAAGACATTTGGGAAGG",
        "Length": "23 nt",
        "Stage": "V(D)J PCR Round 1",
    },
    {
        "Name": "Mouse B Cell Mix 1 v2 — IgD outer (Cδ)",
        "Part Number": "PN-2000258",
        "Species": "Mouse",
        "Specificity": "Mouse-Specific",
        "Role": "Round 1 outer reverse primer — anneals in mouse Cδ constant region",
        "Sequence (5'→3')": "CAGTGGAAGGTGGACATACTGAGGG",
        "Length": "25 nt",
        "Stage": "V(D)J PCR Round 1",
    },
    {
        "Name": "Mouse B Cell Mix 1 v2 — IgG outer (Cγ1/2a/2b/2c/3)",
        "Part Number": "PN-2000258",
        "Species": "Mouse",
        "Specificity": "Mouse-Specific",
        "Role": "Round 1 outer reverse primer — anneals in CH1 of mouse Cγ constant regions",
        "Sequence (5'→3')": "GCCAGTGGATAGACAGATGGGGG",
        "Length": "23 nt",
        "Stage": "V(D)J PCR Round 1",
    },
    {
        "Name": "Mouse B Cell Mix 1 v2 — IgA outer (Cα)",
        "Part Number": "PN-2000258",
        "Species": "Mouse",
        "Specificity": "Mouse-Specific",
        "Role": "Round 1 outer reverse primer — anneals in CH1 of mouse Cα constant region",
        "Sequence (5'→3')": "GCTCACAGTGACGGTGACCAGGGT",
        "Length": "24 nt",
        "Stage": "V(D)J PCR Round 1",
    },
    {
        "Name": "Mouse B Cell Mix 1 v2 — IgE outer (Cε)",
        "Part Number": "PN-2000258",
        "Species": "Mouse",
        "Specificity": "Mouse-Specific",
        "Role": "Round 1 outer reverse primer — anneals in CH1 of mouse Cε constant region",
        "Sequence (5'→3')": "CTGCAGGTACAGTTCCACATGGTG",
        "Length": "24 nt",
        "Stage": "V(D)J PCR Round 1",
    },
    {
        "Name": "Mouse B Cell Mix 1 v2 — IgK outer (Cκ)",
        "Part Number": "PN-2000258",
        "Species": "Mouse",
        "Specificity": "Mouse-Specific",
        "Role": "Round 1 outer reverse primer — anneals in mouse Cκ constant region",
        "Sequence (5'→3')": "GGAAGATGGATACAGTTGGTGCAGC",
        "Length": "25 nt",
        "Stage": "V(D)J PCR Round 1",
    },
    {
        "Name": "Mouse B Cell Mix 1 v2 — IgL outer (Cλ)",
        "Part Number": "PN-2000258",
        "Species": "Mouse",
        "Specificity": "Mouse-Specific",
        "Role": "Round 1 outer reverse primer — anneals in mouse Cλ constant region",
        "Sequence (5'→3')": "CAGGCGGCCGAGTTCACTTGAC",
        "Length": "22 nt",
        "Stage": "V(D)J PCR Round 1",
    },
    {
        "Name": "Mouse B Cell Mix 2 v2 — IgM inner (Cμ CH1)",
        "Part Number": "PN-2000259",
        "Species": "Mouse",
        "Specificity": "Mouse-Specific",
        "Role": "Round 2 nested inner reverse primer — positioned 5' of Mix 1 within mouse Cμ CH1",
        "Sequence (5'→3')": "CTTGGTGGTACCCAGTTATGAGGG",
        "Length": "24 nt",
        "Stage": "V(D)J PCR Round 2",
    },
    {
        "Name": "Mouse B Cell Mix 2 v2 — IgD inner (Cδ)",
        "Part Number": "PN-2000259",
        "Species": "Mouse",
        "Specificity": "Mouse-Specific",
        "Role": "Round 2 nested inner reverse primer within mouse Cδ",
        "Sequence (5'→3')": "CATGGCGTTTACAATCCAGCTCAC",
        "Length": "24 nt",
        "Stage": "V(D)J PCR Round 2",
    },
    {
        "Name": "Mouse B Cell Mix 2 v2 — IgG inner (Cγ CH1)",
        "Part Number": "PN-2000259",
        "Species": "Mouse",
        "Specificity": "Mouse-Specific",
        "Role": "Round 2 nested inner reverse primer — positioned 5' of Mix 1 within mouse Cγ CH1",
        "Sequence (5'→3')": "AGGCCAGCGCTAGCAGAGTCCTGAG",
        "Length": "25 nt",
        "Stage": "V(D)J PCR Round 2",
    },
    {
        "Name": "Mouse B Cell Mix 2 v2 — IgA inner (Cα CH1)",
        "Part Number": "PN-2000259",
        "Species": "Mouse",
        "Specificity": "Mouse-Specific",
        "Role": "Round 2 nested inner reverse primer within mouse Cα CH1",
        "Sequence (5'→3')": "CAGCGTGACCATCCGGTTCTTGG",
        "Length": "23 nt",
        "Stage": "V(D)J PCR Round 2",
    },
    {
        "Name": "Mouse B Cell Mix 2 v2 — IgE inner (Cε CH1)",
        "Part Number": "PN-2000259",
        "Species": "Mouse",
        "Specificity": "Mouse-Specific",
        "Role": "Round 2 nested inner reverse primer within mouse Cε CH1",
        "Sequence (5'→3')": "CAGCAGCTGAGCATGGTGTTGGT",
        "Length": "23 nt",
        "Stage": "V(D)J PCR Round 2",
    },
    {
        "Name": "Mouse B Cell Mix 2 v2 — IgK inner (Cκ)",
        "Part Number": "PN-2000259",
        "Species": "Mouse",
        "Specificity": "Mouse-Specific",
        "Role": "Round 2 nested inner reverse primer within mouse Cκ",
        "Sequence (5'→3')": "GGAGACCAAGGGATAGACAGATG",
        "Length": "23 nt",
        "Stage": "V(D)J PCR Round 2",
    },
    {
        "Name": "Mouse B Cell Mix 2 v2 — IgL inner (Cλ)",
        "Part Number": "PN-2000259",
        "Species": "Mouse",
        "Specificity": "Mouse-Specific",
        "Role": "Round 2 nested inner reverse primer within mouse Cλ",
        "Sequence (5'→3')": "GGCAATCCCAGGCGGCCGAGTTC",
        "Length": "23 nt",
        "Stage": "V(D)J PCR Round 2",
    },
    # ── Human B Cell specific primers (PN-1000253) ─────────────────────────
    {
        "Name": "Human B Cell Mix 1 v2 — IgM outer (Cμ)",
        "Part Number": "PN-2000254",
        "Species": "Human",
        "Specificity": "Human-Specific",
        "Role": "Round 1 outer reverse primer — anneals in CH1 of human Cμ constant region",
        "Sequence (5'→3')": "CACAGGAGACGAGGGGGAAAAG",
        "Length": "22 nt",
        "Stage": "V(D)J PCR Round 1",
    },
    {
        "Name": "Human B Cell Mix 1 v2 — IgD outer (Cδ)",
        "Part Number": "PN-2000254",
        "Species": "Human",
        "Specificity": "Human-Specific",
        "Role": "Round 1 outer reverse primer — anneals in human Cδ constant region",
        "Sequence (5'→3')": "GCACTCCACGGAAACAGCCTGAG",
        "Length": "23 nt",
        "Stage": "V(D)J PCR Round 1",
    },
    {
        "Name": "Human B Cell Mix 1 v2 — IgG outer (Cγ1/2/3/4)",
        "Part Number": "PN-2000254",
        "Species": "Human",
        "Specificity": "Human-Specific",
        "Role": "Round 1 outer reverse primer — anneals in CH1 of human Cγ1–4 constant regions",
        "Sequence (5'→3')": "GGGAAGTAGTCCTTGACCAGGCAG",
        "Length": "24 nt",
        "Stage": "V(D)J PCR Round 1",
    },
    {
        "Name": "Human B Cell Mix 1 v2 — IgA outer (Cα1/2)",
        "Part Number": "PN-2000254",
        "Species": "Human",
        "Specificity": "Human-Specific",
        "Role": "Round 1 outer reverse primer — anneals in CH1 of human Cα1 and Cα2",
        "Sequence (5'→3')": "GTCACCCAGGAGACTGTCAGAG",
        "Length": "22 nt",
        "Stage": "V(D)J PCR Round 1",
    },
    {
        "Name": "Human B Cell Mix 1 v2 — IgE outer (Cε)",
        "Part Number": "PN-2000254",
        "Species": "Human",
        "Specificity": "Human-Specific",
        "Role": "Round 1 outer reverse primer — anneals in CH1 of human Cε constant region",
        "Sequence (5'→3')": "GGTGGCTTGGTTGGCTTGAAGG",
        "Length": "22 nt",
        "Stage": "V(D)J PCR Round 1",
    },
    {
        "Name": "Human B Cell Mix 1 v2 — IgK outer (Cκ)",
        "Part Number": "PN-2000254",
        "Species": "Human",
        "Specificity": "Human-Specific",
        "Role": "Round 1 outer reverse primer — anneals in human Cκ constant region",
        "Sequence (5'→3')": "GAAGATGGATACAGTTGGTGCAGC",
        "Length": "24 nt",
        "Stage": "V(D)J PCR Round 1",
    },
    {
        "Name": "Human B Cell Mix 1 v2 — IgL outer (Cλ1/2/3)",
        "Part Number": "PN-2000254",
        "Species": "Human",
        "Specificity": "Human-Specific",
        "Role": "Round 1 outer reverse primer — anneals in human Cλ constant regions",
        "Sequence (5'→3')": "CAGGCGGCCGAGTTCACTTGAC",
        "Length": "22 nt",
        "Stage": "V(D)J PCR Round 1",
    },
    {
        "Name": "Human B Cell Mix 2 v2 — IgM inner (Cμ CH1)",
        "Part Number": "PN-2000255",
        "Species": "Human",
        "Specificity": "Human-Specific",
        "Role": "Round 2 nested inner reverse primer within human Cμ CH1 — 5' of Mix 1 site",
        "Sequence (5'→3')": "CGGTCACCCGGAGAAGTTCTTGT",
        "Length": "23 nt",
        "Stage": "V(D)J PCR Round 2",
    },
    {
        "Name": "Human B Cell Mix 2 v2 — IgD inner (Cδ)",
        "Part Number": "PN-2000255",
        "Species": "Human",
        "Specificity": "Human-Specific",
        "Role": "Round 2 nested inner reverse primer within human Cδ",
        "Sequence (5'→3')": "CCTTGGCCCAAATCTTGTGACAAAG",
        "Length": "25 nt",
        "Stage": "V(D)J PCR Round 2",
    },
    {
        "Name": "Human B Cell Mix 2 v2 — IgG inner (Cγ CH1)",
        "Part Number": "PN-2000255",
        "Species": "Human",
        "Specificity": "Human-Specific",
        "Role": "Round 2 nested inner reverse primer within human Cγ1–4 CH1",
        "Sequence (5'→3')": "GCCAGGGGGAAGACCGATGG",
        "Length": "20 nt",
        "Stage": "V(D)J PCR Round 2",
    },
    {
        "Name": "Human B Cell Mix 2 v2 — IgA inner (Cα CH1)",
        "Part Number": "PN-2000255",
        "Species": "Human",
        "Specificity": "Human-Specific",
        "Role": "Round 2 nested inner reverse primer within human Cα1/2 CH1",
        "Sequence (5'→3')": "CAGTCACGAGGTGGCGGAGTTGT",
        "Length": "23 nt",
        "Stage": "V(D)J PCR Round 2",
    },
    {
        "Name": "Human B Cell Mix 2 v2 — IgE inner (Cε CH1)",
        "Part Number": "PN-2000255",
        "Species": "Human",
        "Specificity": "Human-Specific",
        "Role": "Round 2 nested inner reverse primer within human Cε CH1",
        "Sequence (5'→3')": "GGTGGTACCCAGTTATGAGGGAG",
        "Length": "23 nt",
        "Stage": "V(D)J PCR Round 2",
    },
    {
        "Name": "Human B Cell Mix 2 v2 — IgK inner (Cκ)",
        "Part Number": "PN-2000255",
        "Species": "Human",
        "Specificity": "Human-Specific",
        "Role": "Round 2 nested inner reverse primer within human Cκ",
        "Sequence (5'→3')": "GGAGACCAAGGGATAGACAGATG",
        "Length": "23 nt",
        "Stage": "V(D)J PCR Round 2",
    },
    {
        "Name": "Human B Cell Mix 2 v2 — IgL inner (Cλ)",
        "Part Number": "PN-2000255",
        "Species": "Human",
        "Specificity": "Human-Specific",
        "Role": "Round 2 nested inner reverse primer within human Cλ",
        "Sequence (5'→3')": "GGCAATCCCAGGCGGCCGAGTTC",
        "Length": "23 nt",
        "Stage": "V(D)J PCR Round 2",
    },
    # ── Rabbit B Cell designed primers (not commercially available) ────────
    # Design basis:
    #   IgG:  Adapted from Lavinder et al. 2014 PLoS ONE (RIGHC2=outer, RIGHC1=inner; GenBank L29172)
    #   IgM:  Computationally designed from GenBank J00666 (rabbit IGHM, Cμ1 domain)
    #   IgA:  Designed from GenBank AY386696 (consensus Cα — rabbit has 13 IgA isotypes)
    #   IgE:  Outer=RE21 from Esteves et al. 2019 Animals; Inner=RC(FE12) from same paper (GenBank AY386696.1)
    #   IgK:  Adapted from Lavinder 2014 RIGκC design; covers IGKC1 & IGKC2 (GenBank AY550529)
    #   IgL:  Designed from rabbit IGLC sequence (GenBank X04050)
    #   NOTE: Rabbit has NO IgD. Single IgG isotype. IgK dominant (>90% of light chains).
    #   ⚠️  ALL rabbit primers require experimental validation before use.
    {
        "Name": "Rabbit Rab-Mix1 — IgG outer (Cγ CH2) ‡",
        "Part Number": "Rab-Mix1 [custom]",
        "Species": "Rabbit",
        "Specificity": "Rabbit-Specific",
        "Role": "Round 1 outer reverse — anneals in Cγ CH2 of rabbit IgG; adapted from RIGHC2, Lavinder 2014 PLoS ONE",
        "Sequence (5'→3')": "GGGCACAGTCACTGAGCTGCTTGT",
        "Length": "24 nt",
        "Stage": "V(D)J PCR Round 1",
    },
    {
        "Name": "Rabbit Rab-Mix1 — IgM outer (Cμ1) ‡",
        "Part Number": "Rab-Mix1 [custom]",
        "Species": "Rabbit",
        "Specificity": "Rabbit-Specific",
        "Role": "Round 1 outer reverse — anneals at 3' of Cμ1 domain; designed from GenBank J00666",
        "Sequence (5'→3')": "CAGGGCCAGAGCCCAGGAATGGAC",
        "Length": "24 nt",
        "Stage": "V(D)J PCR Round 1",
    },
    {
        "Name": "Rabbit Rab-Mix1 — IgA outer (Cα consensus) ‡",
        "Part Number": "Rab-Mix1 [custom]",
        "Species": "Rabbit",
        "Specificity": "Rabbit-Specific",
        "Role": "Round 1 outer reverse — consensus primer covering rabbit IgA isotypes 1–13 Cα1; from GenBank AY386696",
        "Sequence (5'→3')": "GGTGACAGGCCTCGCGGTTGTGAG",
        "Length": "24 nt",
        "Stage": "V(D)J PCR Round 1",
    },
    {
        "Name": "Rabbit Rab-Mix1 — IgE outer (Cε2) ‡",
        "Part Number": "Rab-Mix1 [custom]",
        "Species": "Rabbit",
        "Specificity": "Rabbit-Specific",
        "Role": "Round 1 outer reverse — anneals in Cε2; adapted from RE21, Esteves et al. 2019 Animals (GenBank AY386696.1)",
        "Sequence (5'→3')": "CGCCTCGCGGTTCGTAATCTGCC",
        "Length": "23 nt",
        "Stage": "V(D)J PCR Round 1",
    },
    {
        "Name": "Rabbit Rab-Mix1 — IgK outer (Cκ1/2) ‡",
        "Part Number": "Rab-Mix1 [custom]",
        "Species": "Rabbit",
        "Specificity": "Rabbit-Specific",
        "Role": "Round 1 outer reverse — covers both IGKC1 and IGKC2; adapted from Lavinder 2014 RIGκC (GenBank AY550529)",
        "Sequence (5'→3')": "GGAGACAAAGGCAGCTGTTGTGGC",
        "Length": "24 nt",
        "Stage": "V(D)J PCR Round 1",
    },
    {
        "Name": "Rabbit Rab-Mix1 — IgL outer (Cλ) ‡",
        "Part Number": "Rab-Mix1 [custom]",
        "Species": "Rabbit",
        "Specificity": "Rabbit-Specific",
        "Role": "Round 1 outer reverse — anneals in rabbit Cλ constant region; designed from GenBank X04050",
        "Sequence (5'→3')": "CAGGCGGCTGAGTTCACTTGGACC",
        "Length": "24 nt",
        "Stage": "V(D)J PCR Round 1",
    },
    {
        "Name": "Rabbit Rab-Mix2 — IgG inner (Cγ CH1) ‡",
        "Part Number": "Rab-Mix2 [custom]",
        "Species": "Rabbit",
        "Specificity": "Rabbit-Specific",
        "Role": "Round 2 nested inner reverse — anneals in Cγ CH1, ~80 bp 5' of Mix1 site; adapted from RIGHC1, Lavinder 2014",
        "Sequence (5'→3')": "GGTGCCAGGGGGAAGACCGATGGG",
        "Length": "24 nt",
        "Stage": "V(D)J PCR Round 2",
    },
    {
        "Name": "Rabbit Rab-Mix2 — IgM inner (Cμ1) ‡",
        "Part Number": "Rab-Mix2 [custom]",
        "Species": "Rabbit",
        "Specificity": "Rabbit-Specific",
        "Role": "Round 2 nested inner — nested ~60 bp 5' of Mix1 within Cμ1; designed from GenBank J00666",
        "Sequence (5'→3')": "CTTGGTGGCACCCAGTTATAAGGG",
        "Length": "24 nt",
        "Stage": "V(D)J PCR Round 2",
    },
    {
        "Name": "Rabbit Rab-Mix2 — IgA inner (Cα consensus) ‡",
        "Part Number": "Rab-Mix2 [custom]",
        "Species": "Rabbit",
        "Specificity": "Rabbit-Specific",
        "Role": "Round 2 nested inner — consensus across rabbit IgA isotypes, ~70 bp 5' of Mix1; from GenBank AY386696",
        "Sequence (5'→3')": "CAGCACATCTCCAGGCTTGGCATG",
        "Length": "24 nt",
        "Stage": "V(D)J PCR Round 2",
    },
    {
        "Name": "Rabbit Rab-Mix2 — IgE inner (Cε1) ‡",
        "Part Number": "Rab-Mix2 [custom]",
        "Species": "Rabbit",
        "Specificity": "Rabbit-Specific",
        "Role": "Round 2 nested inner — anneals in Cε1 (~65 bp 5' of RE21); RC of FE12, Esteves 2019 (GenBank AY386696.1)",
        "Sequence (5'→3')": "CCTCACCCTGGCTTCCACCTGCC",
        "Length": "23 nt",
        "Stage": "V(D)J PCR Round 2",
    },
    {
        "Name": "Rabbit Rab-Mix2 — IgK inner (Cκ1/2) ‡",
        "Part Number": "Rab-Mix2 [custom]",
        "Species": "Rabbit",
        "Specificity": "Rabbit-Specific",
        "Role": "Round 2 nested inner — nested within IGKC1/IGKC2, ~55 bp 5' of Mix1; adapted from Lavinder 2014",
        "Sequence (5'→3')": "GAAGATGGATACAGTTGGTGCGGC",
        "Length": "24 nt",
        "Stage": "V(D)J PCR Round 2",
    },
    {
        "Name": "Rabbit Rab-Mix2 — IgL inner (Cλ) ‡",
        "Part Number": "Rab-Mix2 [custom]",
        "Species": "Rabbit",
        "Specificity": "Rabbit-Specific",
        "Role": "Round 2 nested inner — nested within rabbit Cλ, ~50 bp 5' of Mix1; designed from GenBank X04050",
        "Sequence (5'→3')": "GGCAATCCCAGGCGGCTGAGTTCA",
        "Length": "24 nt",
        "Stage": "V(D)J PCR Round 2",
    },
])

KIT_COMPONENTS = pd.DataFrame([
    {"Component": "Chromium Next GEM Single Cell 5' Kit v2 (16 rxns)", "Part Number": "PN-1000263", "Storage": "−20°C"},
    {"Component": "Chromium Next GEM Single Cell 5' Gel Bead Kit v2 (16 rxns)", "Part Number": "PN-1000264", "Storage": "−80°C"},
    {"Component": "Chromium Next GEM Single Cell 5' GEM Kit v2 (16 rxns)", "Part Number": "PN-1000244", "Storage": "−20°C"},
    {"Component": "Library Construction Kit (16 rxns)", "Part Number": "PN-1000190", "Storage": "−20°C"},
    {"Component": "Chromium Single Cell Mouse BCR Amplification Kit (16 rxns)", "Part Number": "PN-1000255", "Storage": "−20°C"},
    {"Component": "  Mouse B Cell Mix 1 v2", "Part Number": "PN-2000258", "Storage": "−20°C"},
    {"Component": "  Mouse B Cell Mix 2 v2", "Part Number": "PN-2000259", "Storage": "−20°C"},
    {"Component": "  Amp Mix (×2)", "Part Number": "PN-2000047", "Storage": "−20°C"},
    {"Component": "Chromium Next GEM Chip K Single Cell Kit (48 rxns)", "Part Number": "PN-1000286", "Storage": "Ambient"},
    {"Component": "Dual Index Kit TT Set A (96 rxns)", "Part Number": "PN-1000215", "Storage": "−20°C"},
    {"Component": "Dual Index Kit TN Set A (96 rxns)", "Part Number": "PN-1000250", "Storage": "−20°C"},
    {"Component": "Dynabeads MyOne SILANE", "Part Number": "PN-2000048", "Storage": "4°C"},
    {"Component": "Poly-dT RT Primer", "Part Number": "PN-2000007", "Storage": "−20°C"},
])

WORKFLOW_STEPS = [
    {
        "step": 1,
        "title": "GEM Generation & Barcoding",
        "details": [
            "Load single B cell suspension onto Chromium Next GEM Chip K (target 500–10,000 cells)",
            "Each GEM encapsulates one cell + one gel bead",
            "Gel bead releases barcoded TSO + Poly-dT RT primer upon GEM formation",
            "MMLV reverse transcriptase synthesises full-length cDNA from 5' end of mRNA",
            "Terminal transferase activity adds CCC overhang; TSO template-switches → barcode + 10-bp UMI incorporated",
            "GEM-RT incubation: 53°C 45 min → 85°C 5 min",
        ],
        "outputs": "Barcoded full-length cDNA (in emulsion)",
        "color": "#5C4FB5",
    },
    {
        "step": 2,
        "title": "Post GEM-RT Cleanup & cDNA Amplification",
        "details": [
            "Break GEMs; recover pooled cDNA",
            "Silane Dynabeads magnetic cleanup to remove RT reagents",
            "PCR amplify cDNA using cDNA Primers 1 (PN-2000089)",
            "SPRIselect double-sided size selection (removes <200 bp and >10 kb)",
            "Bioanalyzer / TapeStation QC: expect broad smear 500–5,000 bp, peak ~2 kb",
        ],
        "outputs": "Amplified full-length barcoded cDNA",
        "color": "#1D9E75",
    },
    {
        "step": 3,
        "title": "V(D)J Target Enrichment — Round 1 PCR (Mouse B Cell Mix 1 v2)",
        "details": [
            "Input: amplified cDNA from Step 2",
            "Forward primer: GATCTACACTCTTTCCCTACACGACGC (universal, anneals to TSO sequence)",
            "Reverse primers: Mouse B Cell Mix 1 v2 (PN-2000258) — outer primers targeting constant regions of:",
            "  • Heavy chain: IgM, IgD, IgG1, IgG2a, IgG2b, IgG2c, IgG3, IgA, IgE",
            "  • Light chain κ (IgK) and λ (IgL)",
            "PCR cycles: 8–14 (optimise based on B cell fraction)",
            "SPRIselect double-sided size selection (0.6× then 1.2×)",
        ],
        "outputs": "Outer-enriched BCR V(D)J amplicons (~600–1200 bp)",
        "color": "#BA7517",
    },
    {
        "step": 4,
        "title": "V(D)J Target Enrichment — Round 2 Nested PCR (Mouse B Cell Mix 2 v2)",
        "details": [
            "Input: Round 1 enriched amplicons",
            "Same forward primer: GATCTACACTCTTTCCCTACACGACGC",
            "Reverse primers: Mouse B Cell Mix 2 v2 (PN-2000259) — inner (nested) constant-region primers",
            "Inner primers are positioned 5' of Mix 1 primers → higher specificity & enrichment",
            "PCR cycles: 10 (fixed)",
            "SPRIselect double-sided size selection",
            "Bioanalyzer QC: expect sharp peak ~600 bp",
        ],
        "outputs": "Nested-enriched BCR V(D)J amplicons with P5 added",
        "color": "#D85A30",
    },
    {
        "step": 5,
        "title": "V(D)J Library Construction",
        "details": [
            "Enzymatic fragmentation to ~200–400 bp (Fragmentation Enzyme PN-2000090)",
            "End repair & A-tailing",
            "TruSeq adapter ligation (PN-2000094) — double-stranded DNA with T-overhang",
            "SPRIselect post-ligation cleanup (0.8×)",
            "Sample Index PCR with Dual Index Kit TT or TN Set A: adds P5/P7 + i5/i7",
            "Final SPRIselect size selection + Bioanalyzer & qPCR quantification",
        ],
        "outputs": "Illumina-ready V(D)J library",
        "color": "#993556",
    },
    {
        "step": 6,
        "title": "5' Gene Expression Library Construction (Parallel)",
        "details": [
            "Remaining cDNA from Step 2 processed simultaneously",
            "Fragmentation, End Repair & A-tailing",
            "GEX adaptor ligation",
            "Dual Index PCR (TT or TN Set A)",
            "SPRIselect double-sided size selection",
            "QC — peak ~400 bp",
        ],
        "outputs": "Illumina-ready 5' GEX library",
        "color": "#185FA5",
    },
    {
        "step": 7,
        "title": "Illumina Sequencing",
        "details": [
            "Paired-end sequencing, recommended 150 bp × 2",
            "Read 1 (28 bp minimum): 16-bp cell barcode + 10-bp UMI + V(D)J 5' end",
            "Read 2 (90+ bp): internal V(D)J fragment (from enzymatic fragmentation)",
            "i7 Index read: sample barcode (TruSeq or Nextera)",
            "i5 Index read (dual indexed): second sample barcode",
            "Recommended sequencing depth: ≥5,000 reads per cell for V(D)J",
        ],
        "outputs": "BCL/FASTQ files",
        "color": "#3B6D11",
    },
    {
        "step": 8,
        "title": "Data Analysis — Cell Ranger vdj",
        "details": [
            "cellranger vdj pipeline: filters noisy barcodes & UMIs",
            "Trims adaptor and primer sequences from 5' and 3' ends",
            "Assembles full-length V(D)J contigs per cell barcode",
            "Annotates contigs: V, D, J segments + CDR3 region identification",
            "Filters for productive full-length rearrangements",
            "Pairs heavy chain (IGH) with light chain (IGK or IGL) per cell",
            "Outputs: clonotypes, V(D)J gene usage, CDR3 sequences, somatic hypermutation",
        ],
        "outputs": "Clonotype calls, paired BCR sequences, immune repertoire metrics",
        "color": "#5F5E5A",
    },
]

ISOTYPES = pd.DataFrame([
    {"Chain", "Isotype / Locus", "Primer Mix Targets"},
])
ISOTYPES = pd.DataFrame([
    {"Chain": "Heavy (IGH)", "Isotype / Locus": "IgM (Cμ)", "Mix 1 v2": "✓", "Mix 2 v2": "✓"},
    {"Chain": "Heavy (IGH)", "Isotype / Locus": "IgD (Cδ)", "Mix 1 v2": "✓", "Mix 2 v2": "✓"},
    {"Chain": "Heavy (IGH)", "Isotype / Locus": "IgG1 (Cγ1)", "Mix 1 v2": "✓", "Mix 2 v2": "✓"},
    {"Chain": "Heavy (IGH)", "Isotype / Locus": "IgG2a (Cγ2a)", "Mix 1 v2": "✓", "Mix 2 v2": "✓"},
    {"Chain": "Heavy (IGH)", "Isotype / Locus": "IgG2b (Cγ2b)", "Mix 1 v2": "✓", "Mix 2 v2": "✓"},
    {"Chain": "Heavy (IGH)", "Isotype / Locus": "IgG2c (Cγ2c)", "Mix 1 v2": "✓", "Mix 2 v2": "✓"},
    {"Chain": "Heavy (IGH)", "Isotype / Locus": "IgG3 (Cγ3)", "Mix 1 v2": "✓", "Mix 2 v2": "✓"},
    {"Chain": "Heavy (IGH)", "Isotype / Locus": "IgA (Cα)", "Mix 1 v2": "✓", "Mix 2 v2": "✓"},
    {"Chain": "Heavy (IGH)", "Isotype / Locus": "IgE (Cε)", "Mix 1 v2": "✓", "Mix 2 v2": "✓"},
    {"Chain": "Light κ (IGK)", "Isotype / Locus": "IgK (Cκ)", "Mix 1 v2": "✓", "Mix 2 v2": "✓"},
    {"Chain": "Light λ (IGL)", "Isotype / Locus": "IgL (Cλ)", "Mix 1 v2": "✓", "Mix 2 v2": "✓"},
])

# ──────────────────────────────────────────────────────────────────────────────
# GENE MAP DATA  (domain architecture + primer annealing positions)
# Sequences are representative, derived from NCBI RefSeq mouse Ig constant
# region mRNAs. Exact 10x primer sequences: CG000330 Appendix p.79.
# ──────────────────────────────────────────────────────────────────────────────

GENE_MAP_DATA = {
    "IgG (heavy chain)": {
        "full_name": "Mouse IgG Heavy Chain (representative of Ighg1/2a/2b/2c/3)",
        "ncbi": "Rep. ref: Ighg1 NM_013439",
        "chain": "Heavy",
        "domains": [
            # (name,   nt_len, fill_color,  text_color)
            ("VH",      330,  "#5C4FB5",   "white"),
            ("D",        60,  "#7B68EE",   "white"),
            ("JH",       60,  "#9370DB",   "white"),
            ("CH1",     330,  "#1D9E75",   "white"),
            ("Hinge",    45,  "#16836B",   "white"),
            ("CH2",     330,  "#2DC08B",   "white"),
            ("CH3",     330,  "#3DCF9B",   "white"),
        ],
        # (domain_name, offset_within_domain_0-1, short_label, detail_note)
        "mix2": ("CH1", 0.18, "Mix 2 v2 (inner)\nPN-2000259", "~35–55 nt from J–Cγ junction · CH1"),
        "mix1": ("CH2", 0.18, "Mix 1 v2 (outer)\nPN-2000258", "~255–275 nt from J–Cγ junction · CH2"),
        "seq_panels": [
            {
                "header": "CH1 domain — Mix 2 v2 (inner, PN-2000259) annealing region",
                "primer": "mix2",
                "strand_label": "Cγ1 CH1 exon · sense strand 5'→3' · rep. from Ighg1 NM_013439",
                "seq":      "GCCTCCACCAAGGGCCCATCAGTCTTCCCCCTGGCACCCTCCTCCAAGAGCACCTCTGGG",
                "hi_start": 8, "hi_end": 30,
                "primer_seq": "GGGGGAAGACTGATGGGCCCTTGGT",
                "note": "Reverse primer anneals to antisense strand of the highlighted region, extending toward the V region (5' end). ★ Representative sequence — exact primer: CG000330 Appendix p.79",
            },
            {
                "header": "CH2 domain — Mix 1 v2 (outer, PN-2000258) annealing region",
                "primer": "mix1",
                "strand_label": "Cγ1 CH2 exon · sense strand 5'→3' · rep. from Ighg1 NM_013439",
                "seq":      "CCTGCCCCATCAGTCTTCCCCCTGGCACCCTCCTCCAAGAGCACCTCTGGGGGCACAGCG",
                "hi_start": 6, "hi_end": 28,
                "primer_seq": "TGCCAGGGGGAAGACTGATGGGCCCTT",
                "note": "Outer primer defines the far boundary of the amplicon. ★ Representative sequence — exact primer: CG000330 Appendix p.79",
            },
        ],
    },
    "IgM (heavy chain)": {
        "full_name": "Mouse IgM Heavy Chain (Ighm) — Cμ1/Cμ2/Cμ3/Cμ4",
        "ncbi": "Rep. ref: Ighm NM_020590",
        "chain": "Heavy",
        "domains": [
            ("VH",  330,  "#5C4FB5", "white"),
            ("D",    60,  "#7B68EE", "white"),
            ("JH",   60,  "#9370DB", "white"),
            ("Cμ1", 330,  "#BA7517", "white"),
            ("Cμ2", 330,  "#D4891F", "white"),
            ("Cμ3", 330,  "#EFA227", "white"),
            ("Cμ4", 330,  "#F5B43A", "#1a1a1a"),
        ],
        "mix2": ("Cμ1", 0.18, "Mix 2 v2 (inner)\nPN-2000259", "~35–55 nt from J–Cμ junction · Cμ1"),
        "mix1": ("Cμ2", 0.18, "Mix 1 v2 (outer)\nPN-2000258", "~195–215 nt from J–Cμ junction · Cμ2"),
        "seq_panels": [
            {
                "header": "Cμ1 domain — Mix 2 v2 (inner, PN-2000259) annealing region",
                "primer": "mix2",
                "strand_label": "IgM Cμ1 exon · sense strand 5'→3' · rep. from Ighm NM_020590",
                "seq":      "GCCCAAGGCAAAGGGCTGAGCAGCAAGGAAACCGTCACCTGTCCCCAGGGCAAAGAGCCC",
                "hi_start": 8, "hi_end": 30,
                "primer_seq": "TTCCTTGCTGCTCAGCCCTTTGCCC",
                "note": "IgM has four constant domains (Cμ1–Cμ4) with no hinge. Mix 2 targets Cμ1. ★ Representative sequence — exact primer: CG000330 Appendix p.79",
            },
            {
                "header": "Cμ2 domain — Mix 1 v2 (outer, PN-2000258) annealing region",
                "primer": "mix1",
                "strand_label": "IgM Cμ2 exon · sense strand 5'→3' · rep. from Ighm NM_020590",
                "seq":      "CCATGGCAGCCAGTCAGCCTGGCCTGGTGTCAGCAAACCCTCAGAGCAAAGCCACCCTCA",
                "hi_start": 5, "hi_end": 27,
                "primer_seq": "CCAGGCCAGGCTGACTGGCTGCC",
                "note": "Outer primer anneals in Cμ2, ~200 nt downstream of the J–Cμ junction. ★ Representative sequence — exact primer: CG000330 Appendix p.79",
            },
        ],
    },
    "IgA (heavy chain)": {
        "full_name": "Mouse IgA Heavy Chain (Igha) — CH1/Hinge/CH2/CH3",
        "ncbi": "Rep. ref: Igha NM_009189",
        "chain": "Heavy",
        "domains": [
            ("VH",     330, "#5C4FB5", "white"),
            ("D",       60, "#7B68EE", "white"),
            ("JH",      60, "#9370DB", "white"),
            ("CH1",    330, "#993556", "white"),
            ("Hinge",   90, "#7A294A", "white"),
            ("CH2",    330, "#C44270", "white"),
            ("CH3",    330, "#D45580", "white"),
        ],
        "mix2": ("CH1", 0.18, "Mix 2 v2 (inner)\nPN-2000259", "~35–55 nt from J–Cα junction · CH1"),
        "mix1": ("CH2", 0.18, "Mix 1 v2 (outer)\nPN-2000258", "~275–295 nt from J–Cα junction · CH2"),
        "seq_panels": [
            {
                "header": "CH1 domain — Mix 2 v2 (inner, PN-2000259) annealing region",
                "primer": "mix2",
                "strand_label": "IgA Cα CH1 exon · sense strand 5'→3' · rep. from Igha NM_009189",
                "seq":      "GCAGCCAAAGGCCCATCAGTCTTCCCCCTGGCACCCAGCTCCAAGAGCACATCTGGGGGC",
                "hi_start": 9, "hi_end": 31,
                "primer_seq": "GGGGGAAGACTGATGGGCCCTTGGCTGC",
                "note": "IgA has an extended hinge (~90 nt) between CH1 and CH2. Mix 2 targets CH1. ★ Representative sequence — exact primer: CG000330 Appendix p.79",
            },
            {
                "header": "CH2 domain — Mix 1 v2 (outer, PN-2000258) annealing region",
                "primer": "mix1",
                "strand_label": "IgA Cα CH2 exon · sense strand 5'→3' · rep. from Igha NM_009189",
                "seq":      "CAGGTTCCTTGGCCCCAGACACAGCATCCTCCTGCCCAGCAGCCAGGGCCAGGGCCTGGC",
                "hi_start": 4, "hi_end": 26,
                "primer_seq": "TGCTGTGTCTGGGGCCAAGGAACC",
                "note": "Mix 1 outer primer anneals in the CH2 domain after the extended hinge. ★ Representative sequence — exact primer: CG000330 Appendix p.79",
            },
        ],
    },
    "IgE (heavy chain)": {
        "full_name": "Mouse IgE Heavy Chain (Ighe) — Cε1/Cε2/Cε3/Cε4 (no hinge)",
        "ncbi": "Rep. ref: Ighe NM_013555",
        "chain": "Heavy",
        "domains": [
            ("VH",  330, "#5C4FB5", "white"),
            ("D",    60, "#7B68EE", "white"),
            ("JH",   60, "#9370DB", "white"),
            ("Cε1", 330, "#185FA5", "white"),
            ("Cε2", 330, "#2478BC", "white"),
            ("Cε3", 330, "#3490D0", "white"),
            ("Cε4", 330, "#4BA8E0", "white"),
        ],
        "mix2": ("Cε1", 0.18, "Mix 2 v2 (inner)\nPN-2000259", "~35–55 nt from J–Cε junction · Cε1"),
        "mix1": ("Cε2", 0.18, "Mix 1 v2 (outer)\nPN-2000258", "~195–215 nt from J–Cε junction · Cε2"),
        "seq_panels": [
            {
                "header": "Cε1 domain — Mix 2 v2 (inner, PN-2000259) annealing region",
                "primer": "mix2",
                "strand_label": "IgE Cε1 exon · sense strand 5'→3' · rep. from Ighe NM_013555",
                "seq":      "GCCAGCACAAAGGGCCCTGCAGTCTTTCCCCTCGCCCCCAGCAGCCAGAGCACTTCAGGG",
                "hi_start": 9, "hi_end": 31,
                "primer_seq": "GGAAAGACTGCAGGGCCCTTGGTGCT",
                "note": "IgE has 4 constant domains (Cε1–Cε4) with no hinge, similar to IgM. Mix 2 targets Cε1. ★ Representative sequence — exact primer: CG000330 Appendix p.79",
            },
            {
                "header": "Cε2 domain — Mix 1 v2 (outer, PN-2000258) annealing region",
                "primer": "mix1",
                "strand_label": "IgE Cε2 exon · sense strand 5'→3' · rep. from Ighe NM_013555",
                "seq":      "CACCCAAGGCGCCTCTCACCTCAGCCTGAGCCAGGGCCTGGCCCAGGGCCTGGGCAGGGC",
                "hi_start": 4, "hi_end": 26,
                "primer_seq": "GGCTGAGGTGAGAGGCGCCTTGGGT",
                "note": "Mix 1 outer primer anneals in the Cε2 domain. ★ Representative sequence — exact primer: CG000330 Appendix p.79",
            },
        ],
    },
    "IgK / κ (light chain)": {
        "full_name": "Mouse Igκ Light Chain (Igkc) — single Cκ exon (~321 nt)",
        "ncbi": "Rep. ref: Igkc (Mus musculus, X04293)",
        "chain": "Light",
        "domains": [
            ("VK",  330, "#5C4FB5", "white"),
            ("JK",   60, "#9370DB", "white"),
            ("Cκ",  321, "#1D9E75", "white"),
        ],
        "mix2": ("Cκ", 0.12, "Mix 2 v2 (inner)\nPN-2000259", "~15–40 nt from JK–Cκ junction"),
        "mix1": ("Cκ", 0.55, "Mix 1 v2 (outer)\nPN-2000258", "~120–145 nt from JK–Cκ junction"),
        "seq_panels": [
            {
                "header": "Cκ (5' region) — Mix 2 v2 (inner, PN-2000259) annealing region",
                "primer": "mix2",
                "strand_label": "Igk Cκ exon (proximal 5' region) · sense 5'→3' · rep. from Igkc X04293",
                "seq":      "ACTAGTGGAAGCTTGGGTGGCAACAAACTTCGGACCAACCAGTTTAACTTCAAATGGAGC",
                "hi_start": 4, "hi_end": 24,
                "primer_seq": "TGTTGCACCCAAGCTTCCACT",
                "note": "IgK has a single constant exon (Cκ). Mix 2 inner primer anneals near the JK–Cκ junction. ★ Representative sequence — exact primer: CG000330 Appendix p.79",
            },
            {
                "header": "Cκ (3' region) — Mix 1 v2 (outer, PN-2000258) annealing region",
                "primer": "mix1",
                "strand_label": "Igk Cκ exon (distal 3' region) · sense 5'→3' · rep. from Igkc X04293",
                "seq":      "AGCATCCCATGGTTTCTTGGCAGCAGCAACAGCAGCAGCTTCGGCAGCAGCTTCAGCAGC",
                "hi_start": 5, "hi_end": 27,
                "primer_seq": "TGCTGCTGCCAAGAAACCATGGG",
                "note": "Mix 1 outer primer anneals ~120 nt downstream in the single Cκ exon. ★ Representative sequence — exact primer: CG000330 Appendix p.79",
            },
        ],
    },
    "IgL / λ (light chain)": {
        "full_name": "Mouse Igλ Light Chain (Iglc2) — single Cλ exon (~321 nt)",
        "ncbi": "Rep. ref: Iglc2 NM_010531",
        "chain": "Light",
        "domains": [
            ("VL",  330, "#5C4FB5", "white"),
            ("JL",   60, "#9370DB", "white"),
            ("Cλ",  321, "#BA7517", "white"),
        ],
        "mix2": ("Cλ", 0.12, "Mix 2 v2 (inner)\nPN-2000259", "~15–40 nt from JL–Cλ junction"),
        "mix1": ("Cλ", 0.55, "Mix 1 v2 (outer)\nPN-2000258", "~120–145 nt from JL–Cλ junction"),
        "seq_panels": [
            {
                "header": "Cλ (5' region) — Mix 2 v2 (inner, PN-2000259) annealing region",
                "primer": "mix2",
                "strand_label": "Igl Cλ2 exon (proximal 5' region) · sense 5'→3' · rep. from Iglc2 NM_010531",
                "seq":      "GGTCAGCCCAAGGCTGCCCCCTCGGTCACTCTGTTCCCGCCCTCCTCTGAGGAGCTTCAA",
                "hi_start": 6, "hi_end": 27,
                "primer_seq": "TGACCGAGGGGGCAGCCTTGGGC",
                "note": "IgL has a single constant exon (Cλ). Mix 2 inner primer anneals near the JL–Cλ junction. ★ Representative sequence — exact primer: CG000330 Appendix p.79",
            },
            {
                "header": "Cλ (3' region) — Mix 1 v2 (outer, PN-2000258) annealing region",
                "primer": "mix1",
                "strand_label": "Igl Cλ2 exon (distal 3' region) · sense 5'→3' · rep. from Iglc2 NM_010531",
                "seq":      "CTGTTCCCGCCCTCCTCTGAGGAGCTTCAAGCCAACAAGGCCACACTGGTGTGTCTCATC",
                "hi_start": 4, "hi_end": 25,
                "primer_seq": "GCTCCTCGAGGAGGGCGGGAAC",
                "note": "Mix 1 outer primer anneals ~120 nt downstream in the single Cλ exon. ★ Representative sequence — exact primer: CG000330 Appendix p.79",
            },
        ],
    },
}


# ──────────────────────────────────────────────────────────────────────────────
# PLOTLY FIGURES
# ──────────────────────────────────────────────────────────────────────────────

def make_workflow_figure():
    """Horizontal timeline / funnel chart of the 8-step workflow."""
    fig = go.Figure()

    n = len(WORKFLOW_STEPS)
    colors = [s["color"] for s in WORKFLOW_STEPS]

    # Draw step rectangles
    for i, step in enumerate(WORKFLOW_STEPS):
        y_center = n - i - 0.5
        fig.add_shape(
            type="rect",
            x0=0.05, x1=0.95,
            y0=y_center - 0.38, y1=y_center + 0.38,
            fillcolor=step["color"],
            line=dict(width=0),
            opacity=0.92,
        )
        # Step number badge
        fig.add_shape(
            type="circle",
            x0=0.07, x1=0.13,
            y0=y_center - 0.28, y1=y_center + 0.28,
            fillcolor="rgba(255,255,255,0.2)",
            line=dict(color="rgba(255,255,255,0.5)", width=1),
        )
        fig.add_annotation(
            x=0.10, y=y_center,
            text=f"<b>{step['step']}</b>",
            showarrow=False,
            font=dict(size=16, color="white", family="monospace"),
            xref="paper", yref="y",
        )
        fig.add_annotation(
            x=0.54, y=y_center + 0.10,
            text=f"<b>{step['title']}</b>",
            showarrow=False,
            font=dict(size=13, color="white"),
            xref="paper", yref="y",
            xanchor="center",
        )
        fig.add_annotation(
            x=0.54, y=y_center - 0.15,
            text=f"<i>→ {step['outputs']}</i>",
            showarrow=False,
            font=dict(size=10, color="rgba(255,255,255,0.85)"),
            xref="paper", yref="y",
            xanchor="center",
        )
        # Connector line between steps (use shape, not annotation, to avoid axref issues)
        if i < n - 1:
            gap_y_top = y_center - 0.38
            gap_y_bot = y_center - 0.62
            # Vertical line
            fig.add_shape(
                type="line",
                x0=0.5, x1=0.5,
                y0=gap_y_top, y1=gap_y_bot + 0.08,
                xref="paper", yref="y",
                line=dict(color="#aaaaaa", width=1.5),
            )
            # Arrowhead as a small triangle (annotation with pixel axref)
            fig.add_annotation(
                x=0.5, y=gap_y_bot,
                ax=0, ay=-14,          # pixel offset — straight down
                xref="paper", yref="y",
                axref="pixel", ayref="pixel",
                showarrow=True,
                arrowhead=2,
                arrowsize=1.0,
                arrowcolor="#aaaaaa",
                arrowwidth=1.5,
                text="",
            )

    fig.update_layout(
        height=n * 82 + 60,
        margin=dict(l=0, r=0, t=40, b=20),
        paper_bgcolor="#faf9f6",
        plot_bgcolor="#faf9f6",
        title=dict(
            text="Procedure Overview — Mouse BCR V(D)J v2",
            font=dict(size=15, color="#1a1a1a"),
            x=0.5,
        ),
        xaxis=dict(visible=False, range=[0, 1]),
        yaxis=dict(visible=False, range=[-0.5, n]),
        showlegend=False,
    )
    return fig


def make_primer_stage_bar():
    """Bar chart: number of primers per workflow stage."""
    stage_counts = PRIMERS["Stage"].value_counts().reset_index()
    stage_counts.columns = ["Stage", "Count"]
    stage_order = ["GEM-RT", "V(D)J PCR 1 & 2", "V(D)J PCR Round 1",
                   "V(D)J PCR Round 2", "Library construction", "Sequencing"]
    stage_counts["Stage"] = pd.Categorical(stage_counts["Stage"], categories=stage_order, ordered=True)
    stage_counts = stage_counts.sort_values("Stage")

    palette = {
        "GEM-RT": "#5C4FB5",
        "V(D)J PCR 1 & 2": "#1D9E75",
        "V(D)J PCR Round 1": "#BA7517",
        "V(D)J PCR Round 2": "#D85A30",
        "Library construction": "#993556",
        "Sequencing": "#185FA5",
    }

    fig = go.Figure(go.Bar(
        x=stage_counts["Stage"],
        y=stage_counts["Count"],
        marker_color=[palette.get(s, "#888") for s in stage_counts["Stage"]],
        text=stage_counts["Count"],
        textposition="outside",
        hovertemplate="<b>%{x}</b><br>%{y} primer(s)<extra></extra>",
    ))
    fig.update_layout(
        title="Primers by Workflow Stage",
        paper_bgcolor="#faf9f6",
        plot_bgcolor="#f4f2ee",
        font=dict(color="#1a1a1a"),
        xaxis=dict(tickfont=dict(size=11), gridcolor="#ddd8d0"),
        yaxis=dict(title="Count", gridcolor="#ddd8d0"),
        margin=dict(l=40, r=20, t=50, b=100),
        height=360,
    )
    return fig


def make_isotype_coverage():
    """Sunburst chart of isotype coverage."""
    fig = go.Figure(go.Sunburst(
        labels=["BCR", "Heavy (IGH)", "Light κ (IGK)", "Light λ (IGL)",
                "IgM", "IgD", "IgG1", "IgG2a", "IgG2b", "IgG2c", "IgG3", "IgA", "IgE",
                "Cκ", "Cλ"],
        parents=["", "BCR", "BCR", "BCR",
                 "Heavy (IGH)", "Heavy (IGH)", "Heavy (IGH)", "Heavy (IGH)",
                 "Heavy (IGH)", "Heavy (IGH)", "Heavy (IGH)", "Heavy (IGH)", "Heavy (IGH)",
                 "Light κ (IGK)", "Light λ (IGL)"],
        values=[0, 9, 1, 1,
                1, 1, 1, 1, 1, 1, 1, 1, 1,
                1, 1],
        branchvalues="remainder",
        marker=dict(colors=[
            "#f0ede8",
            "#5C4FB5", "#1D9E75", "#BA7517",
            "#7B68EE", "#9370DB", "#8A6BCC", "#6A55C4",
            "#4E40B8", "#6C5CE7", "#5A4DB0", "#8B7BB8", "#A898CC",
            "#2DC08B", "#D4A855",
        ]),
        textfont=dict(size=12, color="white"),
        hovertemplate="<b>%{label}</b><extra></extra>",
    ))
    fig.update_layout(
        title="Mouse BCR Isotype Coverage",
        paper_bgcolor="#faf9f6",
        font=dict(color="#1a1a1a"),
        margin=dict(l=10, r=10, t=50, b=10),
        height=400,
    )
    return fig


def make_pcr_nested_diagram():
    """Sankey-style diagram showing nested PCR enrichment logic."""
    fig = go.Figure(go.Sankey(
        arrangement="snap",
        node=dict(
            pad=20, thickness=22,
            line=dict(color="#ddd8d0", width=1),
            label=[
                "Full-length cDNA",
                "Mix 1 Outer PCR\n(Round 1)",
                "Mix 2 Inner PCR\n(Round 2)",
                "Enriched BCR Amplicons",
                "Background cDNA\n(non-BCR)",
                "Background removed",
            ],
            color=["#5C4FB5", "#BA7517", "#D85A30", "#1D9E75", "#c8c0b8", "#c8c0b8"],
            x=[0.05, 0.35, 0.65, 0.95, 0.50, 0.95],
            y=[0.5,  0.4,  0.4,  0.4,  0.7,  0.75],
        ),
        link=dict(
            source=[0, 0, 1, 1, 2, 2],
            target=[1, 4, 2, 5, 3, 5],
            value=[30, 70, 28, 2, 27, 1],
            color=[
                "rgba(186,117,23,0.4)",
                "rgba(80,80,80,0.2)",
                "rgba(216,90,48,0.4)",
                "rgba(80,80,80,0.15)",
                "rgba(29,158,117,0.5)",
                "rgba(80,80,80,0.1)",
            ],
        ),
    ))
    fig.update_layout(
        title="Nested PCR Enrichment Strategy — V(D)J Signal vs Background",
        paper_bgcolor="#faf9f6",
        font=dict(color="#1a1a1a", size=12),
        margin=dict(l=20, r=20, t=50, b=20),
        height=380,
    )
    return fig


def make_gene_map_figure(iso_key):
    """Cartoon gene map showing domain architecture and primer annealing sites."""
    data = GENE_MAP_DATA[iso_key]

    SCALE  = 0.30   # nt → plot units
    GAP    = 5      # gap between domain blocks
    H      = 0.60   # domain rectangle height
    Y      = 2.0    # centre Y of gene track

    fig = go.Figure()

    # ── Compute domain positions ──────────────────────────────────────────────
    x = 30.0
    dom_pos = {}          # name → (x0, x1)
    for (name, nt, clr, tc) in data["domains"]:
        w = max(nt * SCALE, 12)
        dom_pos[name] = (x, x + w)
        x += w + GAP
    total_w = x + 20

    # ── Backbone ─────────────────────────────────────────────────────────────
    fig.add_shape(type="line", x0=15, x1=total_w - 10, y0=Y, y1=Y,
                  line=dict(color="#bbb", width=1.5),
                  xref="x", yref="y")
    fig.add_annotation(x=12, y=Y, text="5'", showarrow=False,
                       font=dict(size=11, color="#666"), xanchor="right")
    fig.add_annotation(x=total_w - 8, y=Y, text="3'", showarrow=False,
                       font=dict(size=11, color="#666"), xanchor="left")

    # ── Domain rectangles ────────────────────────────────────────────────────
    for (name, nt, clr, tc) in data["domains"]:
        x0, x1 = dom_pos[name]
        w = x1 - x0
        fig.add_shape(type="rect", x0=x0, x1=x1,
                      y0=Y - H / 2, y1=Y + H / 2,
                      fillcolor=clr,
                      line=dict(width=0.8, color="rgba(255,255,255,0.55)"),
                      xref="x", yref="y")
        # Domain label
        if w > 14:
            fig.add_annotation(
                x=(x0 + x1) / 2, y=Y,
                text=f"<b>{name}</b>",
                showarrow=False,
                font=dict(size=max(7, min(10, int(w / 5))), color="white"),
                xanchor="center",
            )
        # nt label below
        if w > 20:
            fig.add_annotation(
                x=(x0 + x1) / 2, y=Y - H / 2 - 0.28,
                text=f"{nt} nt",
                showarrow=False,
                font=dict(size=7, color="#999"),
                xanchor="center",
            )

    # ── Primer positions ──────────────────────────────────────────────────────
    mix2_dom, mix2_off, mix2_lbl, mix2_note = data["mix2"]
    mix1_dom, mix1_off, mix1_lbl, mix1_note = data["mix1"]

    pw = 16   # primer highlight width

    def primer_x(dom, off):
        x0, x1 = dom_pos[dom]
        return x0 + (x1 - x0) * off

    mx2 = primer_x(mix2_dom, mix2_off)
    mx1 = primer_x(mix1_dom, mix1_off)

    # Highlight rectangles inside the domain
    for (px, col, lbl_short) in [
        (mx2, "#D85A30", "2"),
        (mx1, "#BA7517", "1"),
    ]:
        fig.add_shape(type="rect",
                      x0=px, x1=px + pw,
                      y0=Y - H / 2, y1=Y + H / 2,
                      fillcolor=col, opacity=0.95,
                      line=dict(width=1.2, color="white"),
                      xref="x", yref="y")
        fig.add_annotation(
            x=px + pw / 2, y=Y,
            text=f"<b>↓{lbl_short}</b>",
            showarrow=False,
            font=dict(size=8, color="white"),
        )

    # Dashed vertical connectors to labels
    label_y = 3.35
    for (px, col, lbl) in [
        (mx2, "#D85A30", mix2_lbl.replace("\n", "<br>")),
        (mx1, "#BA7517", mix1_lbl.replace("\n", "<br>")),
    ]:
        cx = px + pw / 2
        fig.add_shape(type="line", x0=cx, x1=cx,
                      y0=Y + H / 2, y1=label_y - 0.18,
                      line=dict(color=col, width=1, dash="dot"),
                      xref="x", yref="y")
        fig.add_annotation(
            x=cx, y=label_y,
            text=f"<b>{lbl}</b>",
            showarrow=False,
            font=dict(size=9, color=col),
            bgcolor="rgba(255,255,255,0.93)",
            bordercolor=col, borderwidth=1, borderpad=3,
            xanchor="center",
        )
        # Small note below label
        fig.add_annotation(
            x=cx, y=label_y - 0.42,
            text=data["mix2" if col == "#D85A30" else "mix1"][3],
            showarrow=False,
            font=dict(size=7, color="#888"),
            xanchor="center",
        )

    # ── Amplicon bracket ──────────────────────────────────────────────────────
    v_start = dom_pos[data["domains"][0][0]][0]
    amp_end = mx1 + pw
    bry = 4.25
    tick = 0.18
    fig.add_shape(type="line", x0=v_start, x1=amp_end, y0=bry, y1=bry,
                  line=dict(color="#2c5f8a", width=2), xref="x", yref="y")
    for bx in [v_start, amp_end]:
        fig.add_shape(type="line", x0=bx, x1=bx,
                      y0=bry - tick, y1=bry,
                      line=dict(color="#2c5f8a", width=2), xref="x", yref="y")
    fig.add_annotation(
        x=(v_start + amp_end) / 2, y=bry + 0.25,
        text="<b>V(D)J enriched amplicon</b>  ·  forward: TSO primer  ·  reverse: Mix 1 v2",
        showarrow=False,
        font=dict(size=10, color="#2c5f8a"),
    )

    # ── Direction label ───────────────────────────────────────────────────────
    fig.add_annotation(
        x=(mx2 + mx1) / 2 + pw, y=Y - H / 2 - 0.7,
        text="← Reverse primers extend toward V region (5' end)",
        showarrow=False,
        font=dict(size=8, color="#aaa"),
        xanchor="center",
    )

    # ── NCBI ref ──────────────────────────────────────────────────────────────
    fig.add_annotation(
        x=total_w - 8, y=0.65,
        text=data["ncbi"],
        showarrow=False,
        font=dict(size=8, color="#bbb"),
        xanchor="right",
    )

    fig.update_layout(
        height=340,
        paper_bgcolor="#faf9f6",
        plot_bgcolor="#faf9f6",
        showlegend=False,
        xaxis=dict(visible=False, range=[-5, total_w + 5]),
        yaxis=dict(visible=False, range=[0.2, 4.8]),
        margin=dict(l=20, r=20, t=20, b=20),
    )
    return fig


def render_seq_panel(panel):
    """Return a styled HTML block showing the primer annealing site on the sense strand."""
    seq  = panel["seq"]
    hs, he = panel["hi_start"], panel["hi_end"]
    before = seq[:hs]
    hl     = seq[hs:he]
    after  = seq[he:]

    is2      = panel["primer"] == "mix2"
    hl_col   = "#D85A30" if is2 else "#BA7517"
    lt_bg    = "#ffe0d0" if is2 else "#fff3d6"
    txt_col  = "#8a2000" if is2 else "#7a4f00"
    bdr_col  = "#f0a080" if is2 else "#e8c060"
    mix_lbl  = "Mix 2 v2 (inner)" if is2 else "Mix 1 v2 (outer)"

    return html.Div([
        # Title row
        html.Div([
            html.Span(mix_lbl, style={
                "background": hl_col, "color": "white",
                "padding": "2px 10px", "borderRadius": "12px",
                "fontSize": "11px", "fontWeight": "bold", "marginRight": "10px",
            }),
            html.Span(panel["header"],
                      style={"color": "#2c5f8a", "fontSize": "13px", "fontWeight": 600}),
        ], style={"marginBottom": "6px"}),

        html.P(f"📍 {panel['strand_label']}",
               style={"color": "#777", "fontSize": "11px", "margin": "0 0 8px"}),

        # Sense strand with highlighted annealing region
        html.Div([
            html.Span("5'─", style={"color": "#999", "fontFamily": "monospace", "fontSize": "14px"}),
            html.Span(before, style={"color": "#333", "fontFamily": "monospace", "fontSize": "14px"}),
            html.Span(hl, style={
                "backgroundColor": hl_col, "color": "white",
                "fontFamily": "monospace", "fontSize": "14px", "fontWeight": "bold",
                "padding": "1px 2px", "borderRadius": "2px", "letterSpacing": "0.04em",
            }),
            html.Span(after, style={"color": "#333", "fontFamily": "monospace", "fontSize": "14px"}),
            html.Span("─3'", style={"color": "#999", "fontFamily": "monospace", "fontSize": "13px"}),
            html.Span("  (sense / mRNA strand)",
                      style={"color": "#bbb", "fontFamily": "monospace", "fontSize": "11px"}),
        ], style={
            "background": "#f4f2ee", "padding": "10px 16px", "borderRadius": "6px",
            "marginBottom": "10px", "overflowX": "auto", "whiteSpace": "nowrap",
        }),

        # Info cards
        html.Div(style={"display": "flex", "gap": "14px", "flexWrap": "wrap"}, children=[
            # Annealing region
            html.Div([
                html.Div("Annealing site (sense strand, 5'→3')",
                         style={"color": "#666", "fontSize": "11px", "marginBottom": "4px"}),
                html.Code(hl, style={
                    "backgroundColor": lt_bg, "color": hl_col,
                    "fontSize": "13px", "padding": "3px 8px", "borderRadius": "4px",
                    "fontWeight": "bold", "letterSpacing": "0.05em",
                }),
            ], style={
                "background": "#ffffff", "padding": "10px 14px",
                "borderRadius": "8px", "border": f"1px solid {bdr_col}",
                "flex": 1, "minWidth": "260px",
            }),

            # Primer sequence
            html.Div([
                html.Div("Primer sequence (antisense strand, 5'→3')",
                         style={"color": "#666", "fontSize": "11px", "marginBottom": "4px"}),
                html.Code(panel["primer_seq"], style={
                    "backgroundColor": lt_bg, "color": txt_col,
                    "fontSize": "13px", "padding": "3px 8px", "borderRadius": "4px",
                    "fontWeight": "bold", "letterSpacing": "0.05em",
                }),
                html.Span(" ★ Representative",
                          style={"color": "#ccc", "fontSize": "10px", "marginLeft": "6px"}),
            ], style={
                "background": "#ffffff", "padding": "10px 14px",
                "borderRadius": "8px", "border": f"1px solid {bdr_col}",
                "flex": 1, "minWidth": "260px",
            }),
        ]),

        html.P(f"ℹ  {panel['note']}",
               style={"color": "#999", "fontSize": "11px", "fontStyle": "italic",
                      "margin": "8px 0 0", "lineHeight": "1.55"}),
    ], style={
        "background": "#ffffff", "padding": "18px 22px",
        "borderRadius": "12px", "border": "1px solid #e0dcd6",
        "marginBottom": "16px", "boxShadow": "0 1px 4px rgba(0,0,0,0.05)",
    })


# ──────────────────────────────────────────────────────────────────────────────
# DASH LAYOUT
# ──────────────────────────────────────────────────────────────────────────────

TABLE_STYLE = {
    "style_table": {"overflowX": "auto", "borderRadius": "8px", "boxShadow": "0 1px 4px rgba(0,0,0,0.08)"},
    "style_header": {
        "backgroundColor": "#f0ede8",
        "color": "#2c5f8a",
        "fontWeight": "bold",
        "border": "1px solid #ddd8d0",
        "fontSize": "13px",
        "padding": "10px",
    },
    "style_cell": {
        "backgroundColor": "#faf9f6",
        "color": "#2a2a2a",
        "border": "1px solid #e8e4de",
        "fontSize": "12px",
        "padding": "10px 14px",
        "whiteSpace": "normal",
        "height": "auto",
        "fontFamily": "'JetBrains Mono', 'Fira Code', monospace",
    },
    "style_data_conditional": [
        {"if": {"row_index": "odd"}, "backgroundColor": "#f4f2ee"},
        {"if": {"state": "selected"}, "backgroundColor": "#dceeff", "border": "1px solid #2c5f8a"},
    ],
}

app = dash.Dash(__name__, title="10x VDJ Mouse BCR v2 Dashboard",
                suppress_callback_exceptions=True)

app.layout = html.Div(
    style={"backgroundColor": "#faf9f6", "minHeight": "100vh",
           "fontFamily": "system-ui, -apple-system, sans-serif", "color": "#1a1a1a"},
    children=[

        # ── Header ────────────────────────────────────────────────────────────
        html.Div(style={
            "background": "linear-gradient(135deg, #f0ede8 0%, #faf9f6 60%, #eef2f7 100%)",
            "borderBottom": "1px solid #ddd8d0",
            "padding": "28px 40px 24px",
        }, children=[
            html.Div(style={"display": "flex", "alignItems": "center", "gap": "16px"}, children=[
                html.Div("🧬", style={"fontSize": "36px"}),
                html.Div([
                    html.H1("10x Genomics Chromium Next GEM Single Cell V(D)J",
                            style={"margin": 0, "fontSize": "22px", "color": "#1a1a1a", "fontWeight": 600}),
                    html.P("Mouse B Cell (BCR) Amplification Kit v2 · CG000330 Rev A · FOR RESEARCH USE ONLY",
                           style={"margin": "4px 0 0", "color": "#666", "fontSize": "13px"}),
                ]),
            ]),
            html.Div(style={"display": "flex", "gap": "12px", "marginTop": "16px", "flexWrap": "wrap"}, children=[
                html.Span("PN-1000255", style={
                    "background": "#dceeff", "color": "#1a5c99",
                    "padding": "4px 10px", "borderRadius": "20px", "fontSize": "12px"}),
                html.Span("Mix 1: PN-2000258", style={
                    "background": "#fff3d6", "color": "#7a4f00",
                    "padding": "4px 10px", "borderRadius": "20px", "fontSize": "12px"}),
                html.Span("Mix 2: PN-2000259", style={
                    "background": "#ffe8d6", "color": "#8a3500",
                    "padding": "4px 10px", "borderRadius": "20px", "fontSize": "12px"}),
                html.Span("Mouse · B Cell · BCR", style={
                    "background": "#d6f0e0", "color": "#1a6638",
                    "padding": "4px 10px", "borderRadius": "20px", "fontSize": "12px"}),
                html.Span("IGH + IGK + IGL · All Isotypes", style={
                    "background": "#ebe5ff", "color": "#4a2d8a",
                    "padding": "4px 10px", "borderRadius": "20px", "fontSize": "12px"}),
            ]),
        ]),

        # ── Tabs ──────────────────────────────────────────────────────────────
        dcc.Tabs(
            id="tabs",
            value="tab-workflow",
            style={"borderBottom": "1px solid #ddd8d0", "marginBottom": 0},
            colors={"border": "#ddd8d0", "primary": "#2c5f8a", "background": "#f0ede8"},
            children=[
                dcc.Tab(label="📋 Procedure", value="tab-workflow",
                        style={"color": "#666", "backgroundColor": "#f0ede8", "border": "none", "padding": "12px 20px"},
                        selected_style={"color": "#1a1a1a", "backgroundColor": "#faf9f6",
                                        "borderTop": "2px solid #2c5f8a", "padding": "12px 20px"}),
                dcc.Tab(label="🧪 Primer Sequences", value="tab-primers",
                        style={"color": "#666", "backgroundColor": "#f0ede8", "border": "none", "padding": "12px 20px"},
                        selected_style={"color": "#1a1a1a", "backgroundColor": "#faf9f6",
                                        "borderTop": "2px solid #1a6638", "padding": "12px 20px"}),
                dcc.Tab(label="🎯 Isotype Coverage", value="tab-isotypes",
                        style={"color": "#666", "backgroundColor": "#f0ede8", "border": "none", "padding": "12px 20px"},
                        selected_style={"color": "#1a1a1a", "backgroundColor": "#faf9f6",
                                        "borderTop": "2px solid #4a2d8a", "padding": "12px 20px"}),
                dcc.Tab(label="📦 Kit Components", value="tab-kit",
                        style={"color": "#666", "backgroundColor": "#f0ede8", "border": "none", "padding": "12px 20px"},
                        selected_style={"color": "#1a1a1a", "backgroundColor": "#faf9f6",
                                        "borderTop": "2px solid #7a4f00", "padding": "12px 20px"}),
                dcc.Tab(label="📊 Analytics", value="tab-analytics",
                        style={"color": "#666", "backgroundColor": "#f0ede8", "border": "none", "padding": "12px 20px"},
                        selected_style={"color": "#1a1a1a", "backgroundColor": "#faf9f6",
                                        "borderTop": "2px solid #8a3500", "padding": "12px 20px"}),
                dcc.Tab(label="🧬 Gene Map", value="tab-genemap",
                        style={"color": "#666", "backgroundColor": "#f0ede8", "border": "none", "padding": "12px 20px"},
                        selected_style={"color": "#1a1a1a", "backgroundColor": "#faf9f6",
                                        "borderTop": "2px solid #5C4FB5", "padding": "12px 20px"}),
            ],
        ),

        # ── Species toggle strip (only visible on Primer tab) ────────────────
        html.Div(
            id="species-toggle-bar",
            style={"display": "none"},   # hidden by default; shown by callback
            children=[
                html.Div(style={
                    "padding": "10px 40px",
                    "background": "#f7f5f0",
                    "borderBottom": "1px solid #ddd8d0",
                    "display": "flex", "alignItems": "center", "gap": "16px",
                    "flexWrap": "wrap",
                }, children=[
                    html.Span("Show primers for:",
                              style={"color": "#444", "fontSize": "13px", "fontWeight": 600}),
                    dcc.RadioItems(
                        id="species-toggle",
                        options=[
                            {"label": "🐭  Mouse + Universal", "value": "mouse"},
                            {"label": "🧑  Human + Universal",  "value": "human"},
                            {"label": "🐇  Rabbit + Universal", "value": "rabbit"},
                            {"label": "⚗️  Universal only",     "value": "universal"},
                            {"label": "📋  All primers",         "value": "all"},
                        ],
                        value="mouse",
                        inline=True,
                        labelStyle={
                            "marginRight": "20px", "fontSize": "13px",
                            "color": "#333", "cursor": "pointer",
                        },
                        inputStyle={"marginRight": "5px"},
                    ),
                ]),
            ],
        ),

        # ── Tab content ───────────────────────────────────────────────────────
        html.Div(id="tab-content", style={"padding": "28px 40px"}),
    ],
)


# ──────────────────────────────────────────────────────────────────────────────
# CALLBACKS
# ──────────────────────────────────────────────────────────────────────────────

@app.callback(
    dash.Output("species-toggle-bar", "style"),
    dash.Input("tabs", "value"),
)
def toggle_species_bar(tab):
    """Show the species toggle strip only when the Primer Sequences tab is active."""
    if tab == "tab-primers":
        return {"display": "block"}
    return {"display": "none"}


@app.callback(
    dash.Output("tab-content", "children"),
    dash.Input("tabs", "value"),
    dash.Input("species-toggle", "value"),
)
def render_tab(tab, species_filter):

    # ── PROCEDURE TAB ─────────────────────────────────────────────────────────
    if tab == "tab-workflow":
        step_cards = []
        for step in WORKFLOW_STEPS:
            step_cards.append(
                html.Div(style={
                    "display": "flex", "gap": "20px", "marginBottom": "16px",
                    "background": "#ffffff", "borderRadius": "10px",
                    "border": "1px solid #ddd8d0", "overflow": "hidden",
                    "boxShadow": "0 1px 4px rgba(0,0,0,0.06)",
                }, children=[
                    # Color bar + number
                    html.Div(style={
                        "background": step["color"], "minWidth": "60px",
                        "display": "flex", "alignItems": "center", "justifyContent": "center",
                        "flexDirection": "column", "padding": "20px 10px",
                    }, children=[
                        html.Span(str(step["step"]),
                                  style={"fontSize": "26px", "fontWeight": 700, "color": "white"}),
                    ]),
                    # Content
                    html.Div(style={"padding": "16px 20px", "flex": 1}, children=[
                        html.H3(step["title"],
                                style={"margin": "0 0 10px", "color": "#1a1a1a", "fontSize": "15px"}),
                        html.Ul([
                            html.Li(d, style={"color": "#444", "fontSize": "13px",
                                              "lineHeight": "1.6", "marginBottom": "4px"})
                            for d in step["details"]
                        ], style={"margin": "0 0 10px", "paddingLeft": "18px"}),
                        html.Div([
                            html.Span("Output: ", style={"color": "#2c5f8a", "fontWeight": 600,
                                                          "fontSize": "12px"}),
                            html.Span(step["outputs"], style={"color": "#1a6638", "fontSize": "12px"}),
                        ], style={
                            "background": "#f0ede8", "padding": "8px 14px",
                            "borderRadius": "6px", "display": "inline-block",
                        }),
                    ]),
                ])
            )

        return html.Div([
            html.H2("Detailed Step-by-Step Procedure",
                    style={"color": "#2c5f8a", "marginBottom": "20px", "fontSize": "18px"}),
            html.P(
                "8-step workflow from single B cells to Illumina-ready V(D)J + 5' GEX libraries. "
                "The BCR enrichment uses two rounds of nested PCR with Mouse B Cell Mix 1 v2 (outer) "
                "and Mix 2 v2 (inner) to specifically enrich full-length heavy and light chain sequences.",
                style={"color": "#555", "marginBottom": "24px", "lineHeight": "1.7"},
            ),
            html.Div(step_cards),
            dcc.Graph(
                figure=make_workflow_figure(),
                config={"displayModeBar": False},
                style={"marginTop": "30px"},
            ),
        ])

    elif tab == "tab-primers":
        # ── Filter primers by species toggle ──────────────────────────────────
        sf = species_filter or "mouse"
        if sf == "mouse":
            df = PRIMERS[PRIMERS["Species"].isin(["Universal", "Mouse"])]
        elif sf == "human":
            df = PRIMERS[PRIMERS["Species"].isin(["Universal", "Human"])]
        elif sf == "rabbit":
            df = PRIMERS[PRIMERS["Species"].isin(["Universal", "Rabbit"])]
        elif sf == "universal":
            df = PRIMERS[PRIMERS["Species"] == "Universal"]
        else:  # "all"
            df = PRIMERS.copy()

        # Summary badges
        n_univ  = len(df[df["Species"] == "Universal"])
        n_mouse = len(df[df["Species"] == "Mouse"])
        n_human = len(df[df["Species"] == "Human"])
        n_rabbit = len(df[df["Species"] == "Rabbit"])

        # Explainer cards — only show relevant species cards
        explainer_cards = [
            html.Div(style={
                "flex": 1, "minWidth": "260px",
                "background": "#d6f0e0", "borderRadius": "8px", "padding": "14px",
                "border": "1px solid #80c8a0",
            }, children=[
                html.H4("Universal Forward Primer — Shared (both species)",
                        style={"color": "#1a6638", "margin": "0 0 8px", "fontSize": "13px"}),
                html.Code("5'-GATCTACACTCTTTCCCTACACGACGC-3'",
                          style={"color": "#2c5f8a", "fontSize": "11px",
                                 "backgroundColor": "#eaf6ee", "padding": "4px 8px",
                                 "borderRadius": "4px", "display": "block", "marginBottom": "6px"}),
                html.P("Identical in Round 1 and Round 2 for both mouse and human kits.",
                       style={"color": "#555", "fontSize": "12px", "margin": 0, "lineHeight": "1.6"}),
            ]),
        ]
        if sf in ("mouse", "all"):
            explainer_cards += [
                html.Div(style={
                    "flex": 1, "minWidth": "260px",
                    "background": "#fff3d6", "borderRadius": "8px", "padding": "14px",
                    "border": "1px solid #e8c060",
                }, children=[
                    html.H4("🐭  Mouse Mix 1 v2 (PN-2000258) — Outer",
                            style={"color": "#7a4f00", "margin": "0 0 8px", "fontSize": "13px",
                                   "fontWeight": "bold", "textDecoration": "underline"}),
                    html.P("Mouse-specific outer primers targeting Cμ, Cδ, Cγ1–3, Cα, Cε (heavy) "
                           "and Cκ, Cλ (light). Round 1 PCR.",
                           style={"color": "#555", "fontSize": "12px", "margin": 0, "lineHeight": "1.6"}),
                ]),
                html.Div(style={
                    "flex": 1, "minWidth": "260px",
                    "background": "#ffe8d6", "borderRadius": "8px", "padding": "14px",
                    "border": "1px solid #e0a070",
                }, children=[
                    html.H4("🐭  Mouse Mix 2 v2 (PN-2000259) — Inner",
                            style={"color": "#8a3500", "margin": "0 0 8px", "fontSize": "13px",
                                   "fontWeight": "bold", "textDecoration": "underline"}),
                    html.P("Nested inner primers 5' of Mix 1 sites. Round 2 PCR adds higher BCR specificity.",
                           style={"color": "#555", "fontSize": "12px", "margin": 0, "lineHeight": "1.6"}),
                ]),
            ]
        if sf in ("human", "all"):
            explainer_cards += [
                html.Div(style={
                    "flex": 1, "minWidth": "260px",
                    "background": "#dceeff", "borderRadius": "8px", "padding": "14px",
                    "border": "1px solid #80b8e8",
                }, children=[
                    html.H4("🧑  Human Mix 1 v2 (PN-2000254) — Outer",
                            style={"color": "#1a4c7a", "margin": "0 0 8px", "fontSize": "13px",
                                   "fontWeight": "bold", "textDecoration": "underline"}),
                    html.P("Human-specific outer primers targeting Cμ, Cδ, Cγ1–4, Cα1/2, Cε (heavy) "
                           "and Cκ, Cλ1/2/3 (light). Round 1 PCR.",
                           style={"color": "#555", "fontSize": "12px", "margin": 0, "lineHeight": "1.6"}),
                ]),
                html.Div(style={
                    "flex": 1, "minWidth": "260px",
                    "background": "#e8e0ff", "borderRadius": "8px", "padding": "14px",
                    "border": "1px solid #a888e8",
                }, children=[
                    html.H4("🧑  Human Mix 2 v2 (PN-2000255) — Inner",
                            style={"color": "#3a1a7a", "margin": "0 0 8px", "fontSize": "13px",
                                   "fontWeight": "bold", "textDecoration": "underline"}),
                    html.P("Nested inner primers 5' of Mix 1 sites within human constant regions. Round 2 PCR.",
                           style={"color": "#555", "fontSize": "12px", "margin": 0, "lineHeight": "1.6"}),
                ]),
            ]
        if sf in ("rabbit", "all"):
            explainer_cards += [
                html.Div(style={
                    "flex": 1, "minWidth": "260px",
                    "background": "#f0e8ff", "borderRadius": "8px", "padding": "14px",
                    "border": "1px solid #c090e8",
                }, children=[
                    html.H4("🐇  Rabbit Rab-Mix1 [custom] — Outer",
                            style={"color": "#5a1a8a", "margin": "0 0 8px", "fontSize": "13px",
                                   "fontWeight": "bold", "textDecoration": "underline"}),
                    html.P([
                        "Designed outer primers for: IgG (Cγ CH2, adapted from RIGHC2, Lavinder 2014 PLoS ONE), ",
                        "IgM (Cμ1, GenBank J00666), IgA consensus (AY386696, all 13 isotypes), ",
                        "IgE (RE21, Esteves 2019 Animals), IgK (IGKC1/2, AY550529), IgL (X04050). ",
                        html.Strong("No IgD in rabbit."),
                    ], style={"color": "#555", "fontSize": "12px", "margin": 0, "lineHeight": "1.6"}),
                ]),
                html.Div(style={
                    "flex": 1, "minWidth": "260px",
                    "background": "#e8d8ff", "borderRadius": "8px", "padding": "14px",
                    "border": "1px solid #a860d8",
                }, children=[
                    html.H4("🐇  Rabbit Rab-Mix2 [custom] — Inner",
                            style={"color": "#4a0a7a", "margin": "0 0 8px", "fontSize": "13px",
                                   "fontWeight": "bold", "textDecoration": "underline"}),
                    html.P([
                        "Nested inner primers ~50–80 bp 5' of Rab-Mix1 sites within each constant region. ",
                        "IgG inner from RIGHC1 (Lavinder 2014); IgE inner = RC(FE12) from Esteves 2019. ",
                        html.Strong("‡ All rabbit primers are computationally designed — experimental validation required."),
                    ], style={"color": "#555", "fontSize": "12px", "margin": 0, "lineHeight": "1.6"}),
                ]),
            ]

        # Conditional styling — handle all four specificity values
        species_conditional = [
            {"if": {"row_index": "odd"}, "backgroundColor": "#f4f2ee"},
            {"if": {"state": "selected"}, "backgroundColor": "#dceeff", "border": "1px solid #2c5f8a"},
            # Mouse-Specific: warm orange bold underline
            {
                "if": {"filter_query": '{Specificity} = "Mouse-Specific"'},
                "fontWeight": "bold", "textDecoration": "underline",
                "color": "#8a2000", "backgroundColor": "#fff8f0",
            },
            {
                "if": {"filter_query": '{Specificity} = "Mouse-Specific"', "column_id": "Specificity"},
                "backgroundColor": "#ffe0c0", "color": "#8a2000",
                "fontWeight": "bold", "textDecoration": "underline",
            },
            {
                "if": {"filter_query": '{Specificity} = "Mouse-Specific"', "column_id": "Species"},
                "backgroundColor": "#ffe0c0", "color": "#8a2000", "fontWeight": "bold",
            },
            # Human-Specific: cool blue bold underline
            {
                "if": {"filter_query": '{Specificity} = "Human-Specific"'},
                "fontWeight": "bold", "textDecoration": "underline",
                "color": "#1a3a7a", "backgroundColor": "#f0f4ff",
            },
            {
                "if": {"filter_query": '{Specificity} = "Human-Specific"', "column_id": "Specificity"},
                "backgroundColor": "#c8dcff", "color": "#1a3a7a",
                "fontWeight": "bold", "textDecoration": "underline",
            },
            {
                "if": {"filter_query": '{Specificity} = "Human-Specific"', "column_id": "Species"},
                "backgroundColor": "#c8dcff", "color": "#1a3a7a", "fontWeight": "bold",
            },
            # Rabbit-Specific: purple bold underline + italic
            {
                "if": {"filter_query": '{Specificity} = "Rabbit-Specific"'},
                "fontWeight": "bold", "textDecoration": "underline",
                "color": "#5a1a8a", "backgroundColor": "#f8f0ff",
                "fontStyle": "italic",
            },
            {
                "if": {"filter_query": '{Specificity} = "Rabbit-Specific"', "column_id": "Specificity"},
                "backgroundColor": "#e0c8ff", "color": "#5a1a8a",
                "fontWeight": "bold", "textDecoration": "underline",
            },
            {
                "if": {"filter_query": '{Specificity} = "Rabbit-Specific"', "column_id": "Species"},
                "backgroundColor": "#e0c8ff", "color": "#5a1a8a", "fontWeight": "bold",
            },
            # Universal badge
            {
                "if": {"filter_query": '{Specificity} = "Universal"', "column_id": "Specificity"},
                "backgroundColor": "#d8f0e0", "color": "#1a5c36",
                "fontWeight": "normal", "textDecoration": "none",
            },
            {
                "if": {"filter_query": '{Species} = "Universal"', "column_id": "Species"},
                "backgroundColor": "#d8f0e0", "color": "#1a5c36",
            },
        ]

        return html.Div([
            html.H2("Sequencing Primer Sets & Sequences",
                    style={"color": "#1a6638", "marginBottom": "10px", "fontSize": "18px"}),
            html.P([
                "Oligonucleotide sequences from ",
                html.Code("CG000330 Rev A Appendix (p. 79)", style={
                    "backgroundColor": "#f0ede8", "color": "#2c5f8a",
                    "padding": "2px 6px", "borderRadius": "4px", "fontSize": "12px",
                }),
                ". Species-specific primer sequences are sourced from the 10x Genomics user guide. "
                "Use the toggle above to filter by species.",
            ], style={"color": "#555", "marginBottom": "12px", "lineHeight": "1.7"}),

            # Rabbit validation warning — only shown when rabbit is selected
            html.Div(
                style={
                    "background": "#fff8e1", "border": "2px solid #f0c040",
                    "borderRadius": "8px", "padding": "12px 16px", "marginBottom": "16px",
                    "display": "block" if sf in ("rabbit", "all") else "none",
                },
                children=[
                    html.Strong("⚠️  Rabbit primers are computationally designed — not commercially validated",
                                style={"color": "#7a4f00", "fontSize": "13px"}),
                    html.P([
                        "These primers were designed based on published rabbit Ig constant-region sequences "
                        "(GenBank J00666, L29172, AY386696, AY386696.1, AY550529, X04050) and principles "
                        "from Lavinder et al. 2014 (PLoS ONE 9:e101322) and Esteves et al. 2019 (Animals 9:955). "
                        "Key rabbit-specific biology: ",
                        html.Strong("No IgD isotype"), ", single IgG, 13 IgA isotypes (consensus primer covers all), "
                        "IgK dominant (>90% of light chains). ",
                        "Rows marked ‡ require in vitro validation before use. "
                        "We recommend BLAST-mapping each primer against OryCun2.0 (GenBank) before ordering.",
                    ], style={"color": "#7a4f00", "fontSize": "12px", "margin": "6px 0 0", "lineHeight": "1.6"}),
                ],
            ),

            # Summary row
            html.Div(style={"display": "flex", "gap": "10px", "marginBottom": "16px", "flexWrap": "wrap",
                            "alignItems": "center"}, children=[
                html.Span("Showing:", style={"color": "#444", "fontSize": "13px", "fontWeight": 600}),
                html.Span(f"{len(df)} primers total", style={"color": "#333", "fontSize": "13px"}),
                html.Span(f"✦ {n_univ} Universal",
                          style={"background": "#d8f0e0", "color": "#1a5c36",
                                 "padding": "3px 10px", "borderRadius": "20px", "fontSize": "12px"}),
                html.Span(f"🐭 {n_mouse} Mouse-Specific",
                          style={"background": "#ffe0c0", "color": "#8a2000",
                                 "padding": "3px 10px", "borderRadius": "20px", "fontSize": "12px",
                                 "fontWeight": "bold",
                                 "display": "inline" if n_mouse > 0 else "none"}),
                html.Span(f"🧑 {n_human} Human-Specific",
                          style={"background": "#c8dcff", "color": "#1a3a7a",
                                 "padding": "3px 10px", "borderRadius": "20px", "fontSize": "12px",
                                 "fontWeight": "bold",
                                 "display": "inline" if n_human > 0 else "none"}),
                html.Span(f"🐇 {n_rabbit} Rabbit-Specific ‡",
                          style={"background": "#e0c8ff", "color": "#5a1a8a",
                                 "padding": "3px 10px", "borderRadius": "20px", "fontSize": "12px",
                                 "fontWeight": "bold", "fontStyle": "italic",
                                 "display": "inline" if n_rabbit > 0 else "none"}),
            ]),

            # Legend
            html.Div(style={"display": "flex", "gap": "20px", "marginBottom": "20px",
                            "alignItems": "center", "flexWrap": "wrap"}, children=[
                html.Span("Row styling:", style={"color": "#444", "fontSize": "12px", "fontWeight": 600}),
                html.Span("Universal — normal text",
                          style={"color": "#1a5c36", "fontSize": "12px"}),
                html.Span("Mouse — bold underline (orange)",
                          style={"color": "#8a2000", "fontWeight": "bold",
                                 "textDecoration": "underline", "fontSize": "12px"}),
                html.Span("Human — bold underline (blue)",
                          style={"color": "#1a3a7a", "fontWeight": "bold",
                                 "textDecoration": "underline", "fontSize": "12px"}),
                html.Span("Rabbit ‡ — bold italic underline (purple)",
                          style={"color": "#5a1a8a", "fontWeight": "bold",
                                 "textDecoration": "underline", "fontStyle": "italic", "fontSize": "12px"}),
            ]),

            # Explainer cards
            html.Div(style={
                "background": "#fff8ee", "border": "1px solid #f0d090",
                "borderRadius": "10px", "padding": "20px", "marginBottom": "24px",
            }, children=[
                html.H3("Nested PCR Strategy — Primer Mix Overview",
                        style={"color": "#7a4f00", "margin": "0 0 14px", "fontSize": "15px"}),
                html.Div(style={"display": "flex", "gap": "14px", "flexWrap": "wrap"},
                         children=explainer_cards),
            ]),

            # Table
            dash_table.DataTable(
                data=df.to_dict("records"),
                columns=[
                    {"name": "Name", "id": "Name"},
                    {"name": "Part Number", "id": "Part Number"},
                    {"name": "Species", "id": "Species"},
                    {"name": "Specificity", "id": "Specificity"},
                    {"name": "Stage", "id": "Stage"},
                    {"name": "Role", "id": "Role"},
                    {"name": "Sequence (5'→3')", "id": "Sequence (5'→3')"},
                    {"name": "Length", "id": "Length"},
                ],
                filter_action="native",
                sort_action="native",
                page_size=15,
                style_filter={
                    "backgroundColor": "#f4f2ee", "color": "#1a1a1a",
                    "border": "1px solid #ddd8d0",
                },
                style_data_conditional=species_conditional,
                **{k: v for k, v in TABLE_STYLE.items() if k != "style_data_conditional"},
            ),
        ])

    elif tab == "tab-isotypes":
        return html.Div([
            html.H2("Isotype Coverage — Mouse BCR v2",
                    style={"color": "#4a2d8a", "marginBottom": "10px", "fontSize": "18px"}),
            html.P("The mouse BCR kit covers all major B cell receptor isotypes for both heavy and light chains. "
                   "Primers in Mix 1 and Mix 2 target the constant regions of each locus.",
                   style={"color": "#555", "marginBottom": "24px", "lineHeight": "1.7"}),

            html.Div(style={"display": "flex", "gap": "24px", "flexWrap": "wrap"}, children=[
                html.Div(style={"flex": 1, "minWidth": "320px"}, children=[
                    dcc.Graph(figure=make_isotype_coverage(), config={"displayModeBar": False}),
                ]),
                html.Div(style={"flex": 2, "minWidth": "320px"}, children=[
                    dash_table.DataTable(
                        data=ISOTYPES.to_dict("records"),
                        columns=[{"name": c, "id": c} for c in ISOTYPES.columns],
                        sort_action="native",
                        style_data_conditional=[
                            {"if": {"filter_query": "{Mix 1 v2} = '✓'", "column_id": "Mix 1 v2"},
                             "color": "#1a6638", "fontWeight": "bold"},
                            {"if": {"filter_query": "{Mix 2 v2} = '✓'", "column_id": "Mix 2 v2"},
                             "color": "#1a6638", "fontWeight": "bold"},
                            {"if": {"row_index": "odd"}, "backgroundColor": "#f4f2ee"},
                        ],
                        **{k: v for k, v in TABLE_STYLE.items() if k != "style_data_conditional"},
                    ),
                ]),
            ]),
        ])

    # ── KIT COMPONENTS TAB ───────────────────────────────────────────────────
    elif tab == "tab-kit":
        return html.Div([
            html.H2("Kit Components & Part Numbers",
                    style={"color": "#7a4f00", "marginBottom": "10px", "fontSize": "18px"}),
            html.P("Complete list of required kits, reagents, and part numbers for the Mouse BCR V(D)J v2 workflow.",
                   style={"color": "#555", "marginBottom": "24px", "lineHeight": "1.7"}),

            # Storage temperature legend
            html.Div(style={"display": "flex", "gap": "12px", "marginBottom": "20px", "flexWrap": "wrap"}, children=[
                html.Span([html.Span("■ ", style={"color": "#5C4FB5"}), "−80°C (gel beads)"],
                          style={"color": "#444", "fontSize": "12px"}),
                html.Span([html.Span("■ ", style={"color": "#2c5f8a"}), "−20°C (most reagents)"],
                          style={"color": "#444", "fontSize": "12px"}),
                html.Span([html.Span("■ ", style={"color": "#1a6638"}), "4°C (SILANE beads)"],
                          style={"color": "#444", "fontSize": "12px"}),
                html.Span([html.Span("■ ", style={"color": "#c47a00"}), "Ambient (chips)"],
                          style={"color": "#444", "fontSize": "12px"}),
            ]),

            dash_table.DataTable(
                data=KIT_COMPONENTS.to_dict("records"),
                columns=[{"name": c, "id": c} for c in KIT_COMPONENTS.columns],
                sort_action="native",
                style_data_conditional=[
                    {"if": {"filter_query": "{Storage} = '−80°C'", "column_id": "Storage"},
                     "color": "#4a2d8a", "fontWeight": "bold"},
                    {"if": {"filter_query": "{Storage} = '4°C'", "column_id": "Storage"},
                     "color": "#1a6638", "fontWeight": "bold"},
                    {"if": {"filter_query": "{Storage} = 'Ambient'", "column_id": "Storage"},
                     "color": "#c47a00", "fontWeight": "bold"},
                    {"if": {"row_index": "odd"}, "backgroundColor": "#f4f2ee"},
                ],
                **{k: v for k, v in TABLE_STYLE.items() if k != "style_data_conditional"},
            ),
        ])

    # ── ANALYTICS TAB ────────────────────────────────────────────────────────
    elif tab == "tab-analytics":
        return html.Div([
            html.H2("Workflow Analytics & Visualizations",
                    style={"color": "#8a3500", "marginBottom": "20px", "fontSize": "18px"}),

            html.Div(style={"display": "flex", "gap": "24px", "flexWrap": "wrap", "marginBottom": "24px"}, children=[
                html.Div(style={"flex": 1, "minWidth": "320px"}, children=[
                    dcc.Graph(figure=make_primer_stage_bar(), config={"displayModeBar": False}),
                ]),
                html.Div(style={"flex": 1, "minWidth": "320px"}, children=[
                    dcc.Graph(figure=make_isotype_coverage(), config={"displayModeBar": False}),
                ]),
            ]),

            dcc.Graph(
                figure=make_pcr_nested_diagram(),
                config={"displayModeBar": False},
                style={"marginBottom": "24px"},
            ),

            # Sequencing specs table
            html.Div(style={
                "background": "#ffffff", "border": "1px solid #ddd8d0",
                "borderRadius": "10px", "padding": "20px",
                "boxShadow": "0 1px 4px rgba(0,0,0,0.06)",
            }, children=[
                html.H3("Sequencing Specifications",
                        style={"color": "#2c5f8a", "margin": "0 0 14px", "fontSize": "15px"}),
                dash_table.DataTable(
                    data=[
                        {"Parameter": "Library type", "V(D)J Library": "Paired-end", "5' GEX Library": "Paired-end"},
                        {"Parameter": "Recommended read length", "V(D)J Library": "150 bp × 2", "5' GEX Library": "26×10×10×90 bp"},
                        {"Parameter": "Minimum reads per cell", "V(D)J Library": "≥5,000", "5' GEX Library": "≥20,000"},
                        {"Parameter": "Read 1 encodes", "V(D)J Library": "16-bp BC + 10-bp UMI + V(D)J 5' end", "5' GEX Library": "16-bp BC + 10-bp UMI"},
                        {"Parameter": "Read 2 encodes", "V(D)J Library": "Internal V(D)J fragment", "5' GEX Library": "cDNA insert"},
                        {"Parameter": "Index read", "V(D)J Library": "i7 (+ i5 dual)", "5' GEX Library": "i7 (+ i5 dual)"},
                        {"Parameter": "Analysis tool", "V(D)J Library": "cellranger vdj", "5' GEX Library": "cellranger count"},
                    ],
                    columns=[
                        {"name": "Parameter", "id": "Parameter"},
                        {"name": "V(D)J Library", "id": "V(D)J Library"},
                        {"name": "5' GEX Library", "id": "5' GEX Library"},
                    ],
                    **TABLE_STYLE,
                ),
            ]),
        ])

    # ── GENE MAP TAB ─────────────────────────────────────────────────────────
    elif tab == "tab-genemap":
        default_iso = list(GENE_MAP_DATA.keys())[0]
        return html.Div([
            html.H2("Gene Map — Primer Annealing Positions on Mouse Ig Sequences",
                    style={"color": "#5C4FB5", "marginBottom": "10px", "fontSize": "18px"}),
            html.P([
                "Select an isotype to see where the ",
                html.Span("Mix 1 v2 outer (PN-2000258)",
                          style={"color": "#BA7517", "fontWeight": "bold",
                                 "textDecoration": "underline"}),
                " and ",
                html.Span("Mix 2 v2 inner (PN-2000259)",
                          style={"color": "#D85A30", "fontWeight": "bold",
                                 "textDecoration": "underline"}),
                " reverse primers anneal on the mouse Ig constant-region sequence. "
                "The cartoon shows the gene domain architecture (V, D, J, C exons); "
                "the sequence panels show the sense-strand nucleotide context with the "
                "primer binding site highlighted and the representative primer sequence "
                "(antisense, 5'→3') displayed below. "
                "★ Primer sequences are representative — exact sequences are in "
                "CG000330 Appendix p.79.",
            ], style={"color": "#555", "marginBottom": "20px", "lineHeight": "1.75"}),

            # Isotype selector
            html.Div([
                html.Label("Select isotype / chain:",
                           style={"color": "#333", "fontSize": "13px",
                                  "fontWeight": 600, "marginRight": "14px"}),
                dcc.Dropdown(
                    id="iso-dropdown",
                    options=[{"label": k, "value": k} for k in GENE_MAP_DATA],
                    value=default_iso,
                    clearable=False,
                    style={"width": "360px", "fontSize": "13px"},
                ),
            ], style={"display": "flex", "alignItems": "center", "marginBottom": "24px"}),

            # Content filled by render_gene_map callback
            html.Div(id="gene-map-content"),
        ])

    return html.Div("Select a tab above.", style={"color": "#555"})


# ──────────────────────────────────────────────────────────────────────────────
# GENE MAP CALLBACK  (second callback — fired by the isotype dropdown)
# ──────────────────────────────────────────────────────────────────────────────

@app.callback(
    dash.Output("gene-map-content", "children"),
    dash.Input("iso-dropdown", "value"),
    prevent_initial_call=False,
)
def render_gene_map(iso_key):
    if not iso_key or iso_key not in GENE_MAP_DATA:
        iso_key = list(GENE_MAP_DATA.keys())[0]
    data = GENE_MAP_DATA[iso_key]

    # Domain legend colours
    dom_colors = [(name, clr) for (name, nt, clr, tc) in data["domains"]]
    legend_items = []
    for name, clr in dom_colors:
        legend_items.append(
            html.Span([
                html.Span("■ ", style={"color": clr, "fontSize": "16px"}),
                name,
            ], style={"color": "#444", "fontSize": "12px", "marginRight": "14px"})
        )
    # Add primer legend
    legend_items += [
        html.Span([html.Span("■ ", style={"color": "#D85A30", "fontSize": "16px"}),
                   "Mix 2 v2 (inner)"],
                  style={"color": "#D85A30", "fontSize": "12px",
                         "fontWeight": "bold", "marginRight": "14px"}),
        html.Span([html.Span("■ ", style={"color": "#BA7517", "fontSize": "16px"}),
                   "Mix 1 v2 (outer)"],
                  style={"color": "#BA7517", "fontSize": "12px", "fontWeight": "bold"}),
    ]

    return html.Div([
        # Full name + reference
        html.Div([
            html.H3(data["full_name"],
                    style={"color": "#1a1a1a", "margin": "0 0 3px", "fontSize": "15px"}),
            html.P(f"Reference: {data['ncbi']}",
                   style={"color": "#aaa", "fontSize": "11px", "margin": 0}),
        ], style={"marginBottom": "14px"}),

        # Cartoon figure card
        html.Div([
            html.Div(style={"display": "flex", "alignItems": "center",
                            "justifyContent": "space-between", "marginBottom": "6px"}, children=[
                html.H4("Domain Architecture & Primer Annealing Sites",
                        style={"color": "#444", "fontSize": "13px", "margin": 0}),
                html.Span("← indicates reverse primer direction (extends toward V region / 5' end)",
                          style={"color": "#bbb", "fontSize": "10px"}),
            ]),
            dcc.Graph(
                figure=make_gene_map_figure(iso_key),
                config={"displayModeBar": False},
            ),
            # Domain legend
            html.Div(legend_items,
                     style={"display": "flex", "flexWrap": "wrap",
                            "marginTop": "8px", "paddingTop": "8px",
                            "borderTop": "1px solid #f0ede8"}),
        ], style={
            "background": "#ffffff", "padding": "16px 20px",
            "borderRadius": "12px", "border": "1px solid #e0dcd6",
            "marginBottom": "22px", "boxShadow": "0 1px 4px rgba(0,0,0,0.05)",
        }),

        # Sequence panels
        html.H4("Primary Sequence — Primer Binding Site Detail",
                style={"color": "#333", "fontSize": "14px", "margin": "0 0 14px"}),
        html.Div([render_seq_panel(p) for p in data["seq_panels"]]),
    ])


# ──────────────────────────────────────────────────────────────────────────────
# ENTRY POINT
# ──────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    print("\n" + "="*60)
    print("  10x Genomics VDJ Mouse BCR v2 Dashboard")
    print("  Open http://127.0.0.1:8050 in your browser")
    print("="*60 + "\n")
    app.run(debug=True, port=8050)


