#!/usr/bin/env python3
"""find_guides.py — Scan ABL1 early coding exons for SpCas9 guide RNAs.

SpCas9 guide = 20 nt protospacer immediately 5' of an NGG PAM.
We scan BOTH strands of the target region (exons 2-3, the early coding exons).
"""
import os, re
from Bio import SeqIO
from Bio.Seq import Seq

HERE = os.path.dirname(os.path.abspath(__file__))
rec = SeqIO.read(os.path.join(HERE,"data","ABL1_NM_005157.gb"),"genbank")
mrna = str(rec.seq)

# Target window: exon 2 + exon 3 (early, fully-coding) = mRNA 272..742
TARGET_START, TARGET_END = 272, 742
target = mrna[TARGET_START:TARGET_END]
print(f"Target region: mRNA {TARGET_START}-{TARGET_END} ({len(target)} bp, ABL1 exons 2-3)\n")

def find_guides_on_strand(seq, strand, offset):
    """Find all 20nt + NGG (PAM) guides on the given strand."""
    guides = []
    # NGG PAM: look for GG, guide is the 20nt upstream
    for m in re.finditer(r'(?=([ACGT]{20})([ACGT]GG))', seq):
        protospacer = m.group(1)
        pam = m.group(2)
        pos = m.start()
        guides.append({
            "guide": protospacer, "pam": pam, "strand": strand,
            "mrna_pos": offset + pos if strand == "+" else offset + (len(seq) - pos - 23),
            "gc": round(100*(protospacer.count("G")+protospacer.count("C"))/20, 1),
        })
    return guides

# Plus strand
plus = find_guides_on_strand(target, "+", TARGET_START)
# Minus strand (reverse complement)
target_rc = str(Seq(target).reverse_complement())
minus = find_guides_on_strand(target_rc, "-", TARGET_START)

all_guides = plus + minus
print(f"Candidate guides found: {len(plus)} (+ strand) + {len(minus)} (- strand) = {len(all_guides)}\n")

# Basic quality filters every CRISPR designer applies:
# - GC content 40-70% (too low = weak binding, too high = off-targets)
# - no poly-T (TTTT terminates Pol III transcription of the guide)
# - no extreme homopolymers
def passes(g):
    if not (40 <= g["gc"] <= 75): return False
    if "TTTT" in g["guide"]: return False
    return True

filtered = [g for g in all_guides if passes(g)]
print(f"After quality filters (GC 40-75%, no poly-T): {len(filtered)} guides\n")

print("Sample of candidate guides:")
print(f"{'guide (20nt)':22s} {'PAM':4s} {'strand':6s} {'GC%':5s} {'pos':5s}")
for g in filtered[:10]:
    print(f"{g['guide']:22s} {g['pam']:4s} {g['strand']:6s} {g['gc']:5} {g['mrna_pos']}")

import json
with open(os.path.join(HERE,"data","candidate_guides.json"),"w") as f:
    json.dump(filtered, f, indent=2)
print(f"\nSaved {len(filtered)} candidate guides to data/candidate_guides.json")
