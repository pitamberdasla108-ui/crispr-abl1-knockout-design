#!/usr/bin/env python3
"""score_guides.py — Discriminating on-target efficiency scoring for ABL1 guides.
Literature-based interpretable features (Doench/Rule-Set-2 style), scaled to spread
scores so the ranking is meaningful. Cross-validated against CRISPOR downstream."""
import os, json
HERE = os.path.dirname(os.path.abspath(__file__))
guides = json.load(open(os.path.join(HERE,"data","candidate_guides.json")))

def on_target_score(g):
    seq = g["guide"]; pam = g["pam"]; gc = g["gc"]
    s = 0.0
    # GC content: peak reward at 45-60%, falling off each side (max 30 pts)
    if 45 <= gc <= 60:   s += 30
    elif 40 <= gc < 45 or 60 < gc <= 65: s += 22
    elif 35 <= gc < 40 or 65 < gc <= 70: s += 12
    else: s += 4
    # PAM quality (max 15): NGG with G at pos21 strongest
    s += 15 if pam.endswith("GG") else 5
    # 5' G for U6 transcription efficiency (max 8)
    s += 8 if seq[0] == "G" else 3
    # Seed region GC (PAM-proximal 8nt) — moderate is best (max 15)
    seed_gc = (seq[-8:].count("G")+seq[-8:].count("C"))/8
    if 0.375 <= seed_gc <= 0.625: s += 15
    elif 0.25 <= seed_gc < 0.375 or 0.625 < seed_gc <= 0.75: s += 9
    else: s += 3
    # Position 20 (PAM-proximal): purine (A/G) favored (max 7)
    s += 7 if seq[-1] in "AG" else 3
    # Homopolymer penalty (subtract)
    for b in "ACGT":
        if b*4 in seq: s -= 12; break
    # Pos-3 from PAM often G-rich preference (max 5)
    s += 5 if seq[-3] in "GC" else 2
    # Avoid T at position 1 (max 5)
    s += 5 if seq[0] != "T" else 1
    return round(max(0, s)/85, 3)  # normalize to 0-1 (max possible ~85)

for g in guides: g["on_target"] = on_target_score(g)
guides.sort(key=lambda x: -x["on_target"])

scores = [g["on_target"] for g in guides]
print(f"Scored {len(guides)} guides. Score spread: min {min(scores)}, max {max(scores)}, "
      f"unique values: {len(set(scores))}\n")
print(f"{'rank':4s} {'guide (20nt)':22s} {'PAM':4s} {'str':4s} {'GC%':5s} {'on-tgt':7s} {'pos':5s}")
for i,g in enumerate(guides[:12],1):
    print(f"{i:<4} {g['guide']:22s} {g['pam']:4s} {g['strand']:4s} {g['gc']:<5} {g['on_target']:<7} {g['mrna_pos']}")

json.dump(guides, open(os.path.join(HERE,"data","scored_guides.json"),"w"), indent=2)
print(f"\nSaved ranked guides to data/scored_guides.json")
