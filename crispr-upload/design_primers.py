#!/usr/bin/env python3
"""design_primers.py — Design PCR primers flanking the CRISPR cut site to
amplify the edited ABL1 region for knockout validation (T7E1 / Sanger / TIDE).
Note: ABL1's 5' coding region is GC-rich, so GC window allows up to 65%."""
import os, json
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction, MeltingTemp as mt
from Bio.Seq import Seq

HERE = os.path.dirname(os.path.abspath(__file__))
rec = SeqIO.read(os.path.join(HERE,"data","ABL1_NM_005157.gb"),"genbank")
mrna = str(rec.seq)

GUIDE = "GTCTGAGTGAAGCCGCTCGT"          # CRISPOR top pick, MIT 98
gpos = mrna.find(GUIDE)
cut_site = gpos + 17                      # Cas9 cuts 3bp 5' of PAM
print(f"Selected guide: {GUIDE}  (CRISPOR MIT specificity 98)")
print(f"Guide at mRNA {gpos}, predicted cut site at mRNA {cut_site}\n")

def tm(s): return round(mt.Tm_NN(Seq(s)),1)
def gc(s): return round(gc_fraction(s)*100,1)

def pick(seq, center, direction, win=40):
    best=None
    for start in range(center-win, center+win):
        for L in range(18,26):
            if direction=="fw":
                p = seq[start:start+L]; spos=start
            else:
                p = str(Seq(seq[start:start+L]).reverse_complement()); spos=start+L
            if len(p)<18: continue
            t=tm(p); g=gc(p)
            if 56<=t<=64 and 40<=g<=65:
                score = abs(t-60)+abs(g-52)*0.05
                if best is None or score<best[0]:
                    best=(score,p,t,g,spos)
    return best

fw = pick(mrna, cut_site-160, "fw")
rv = pick(mrna, cut_site+160, "rev")

print("=== PCR primers for knockout validation ===\n")
_,fp,ft,fg,fs = fw
print(f"Forward: 5'-{fp}-3'   (len {len(fp)}, Tm {ft}C, GC {fg}%, pos {fs})")
_,rp,rt,rg,rs = rv
print(f"Reverse: 5'-{rp}-3'   (len {len(rp)}, Tm {rt}C, GC {rg}%, pos {rs})")

amp_start, amp_end = fs, rs
amplicon = mrna[amp_start:amp_end]
left = cut_site-amp_start; right = amp_end-cut_site
print(f"\nAmplicon: {len(amplicon)} bp")
print(f"Cut site sits {left} bp from 5' end / {right} bp from 3' end")
print(f"T7E1 assay -> two fragments of ~{left} bp and ~{right} bp on the gel")
print(f"(unedited DNA stays one {len(amplicon)} bp band; edited shows cleavage)")

json.dump({"guide":GUIDE,"cut_site":cut_site,"fw_primer":fp,"rv_primer":rp,
           "amplicon_bp":len(amplicon),"t7e1_fragments":[left,right]},
          open(os.path.join(HERE,"data","validation_primers.json"),"w"), indent=2)
print("\nSaved data/validation_primers.json")
print("\nNOTE: designed on mRNA/CDS coordinates for this demonstration. A real")
print("wet-lab assay requires primers on GENOMIC DNA (introns included) and a")
print("BLAT specificity check — documented as the next step.")
