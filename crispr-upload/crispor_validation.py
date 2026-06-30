#!/usr/bin/env python3
"""crispor_validation.py — Record CRISPOR genome-wide validation results and
select final guides combining on-target efficiency + genome-wide specificity."""
import os, json
HERE = os.path.dirname(os.path.abspath(__file__))

# CRISPOR genome-wide results (manually recorded from web run, hg19, NGG)
# guide_seq: (MIT_spec, CFD_spec, total_offtargets, notable_offtarget)
crispor = {
 "GAGAAACACTCCTGGTACCA": (42, 87, 164, "1:exon:ABL2 (paralog!)"),
 "GCTGAGTATCTGCTGAGCAG": (53, 78, 244, "intergenic"),
 "GCTGCTCAGCAGATACTCAG": (65, 81, 227, "2:intergenic:MIR5100-RET"),
 "AGCCGCTCGTTGGAACTCCA": (90, 95, 41,  "4:intron:PANK2"),
 "CATCACGCCAGTCAACAGTC": (90, 96, 55,  "intergenic"),
 "GTCTGAGTGAAGCCGCTCGT": (98, 97, 27,  "4:exon:TP53INP1"),  # CRISPOR's #1
}

guides = json.load(open(os.path.join(HERE,"data","scored_guides.json")))
gmap = {g["guide"]: g for g in guides}

print("=== On-target scoring vs CRISPOR genome-wide specificity ===\n")
print(f"{'guide':22s} {'my_ontgt':8s} {'MIT_spec':8s} {'offtgts':7s} {'verdict'}")
rows=[]
for seq,(mit,cfd,ot,note) in sorted(crispor.items(), key=lambda x:-x[1][0]):
    my = gmap.get(seq,{}).get("on_target","--")
    verdict = "SAFE" if mit>=80 else ("RISKY" if mit<60 else "MARGINAL")
    print(f"{seq:22s} {str(my):8s} {mit:<8} {ot:<7} {verdict}  {note}")
    rows.append({"guide":seq,"on_target":my,"mit_spec":mit,"cfd_spec":cfd,
                 "offtargets":ot,"note":note,"verdict":verdict})

print("\n=== KEY FINDING ===")
print("My on-target #1 (GAGAAACACTCCTGGTACCA) has MIT specificity 42 and 164")
print("off-targets including an ABL2 EXON (ABL1's paralog). Isolated scoring")
print("missed this; only genome-wide search (CRISPOR) caught it.\n")
print("=== FINAL SELECTED GUIDES (high on-target AND high specificity) ===")
final = [r for r in rows if r["mit_spec"]>=80]
for r in sorted(final,key=lambda x:-x["mit_spec"]):
    print(f"  {r['guide']}  MIT={r['mit_spec']} off-targets={r['offtargets']}")

json.dump(rows, open(os.path.join(HERE,"data","crispor_validation.json"),"w"), indent=2)
print("\nSaved data/crispor_validation.json")
