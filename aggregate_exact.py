"""Aggregate the exact stage's cost and accuracy out of the parity CSV.

Run: python aggregate_exact.py  (from the repo root)
"""

import csv
import statistics as st
from collections import defaultdict

PATH = "target/parity_geos_report.csv"


def num(v):
    if v is None or v == "":
        return None
    try:
        return float(v)
    except ValueError:
        return None


rows = list(csv.DictReader(open(PATH, encoding="utf-8")))
print("cases in csv:", len(rows))

by_kind = defaultdict(list)
for r in rows:
    by_kind[r["kind"]].append(r)

def pct(vals, p):
    if not vals:
        return float("nan")
    v = sorted(vals)
    i = min(len(v) - 1, max(0, int(round(p / 100.0 * (len(v) - 1)))))
    return v[i]

print("\n=== exact stage, per case, by corpus kind ===")
print(f"{'kind':<10} {'n':>5} {'total s':>9} {'median us':>10} {'p90 us':>10} {'p99 us':>10} {'max us':>11} {'max cand':>9}")
for kind, rs in sorted(by_kind.items()):
    us = [num(r["us_exact"]) for r in rs]
    us = [u for u in us if u is not None]
    cands = [num(r["candidates_exact"]) or 0 for r in rs]
    print(
        f"{kind:<10} {len(rs):>5} {sum(us) / 1e6:>9.3f} {st.median(us):>10.1f} "
        f"{pct(us, 90):>10.1f} {pct(us, 99):>10.1f} {max(us):>11.1f} {max(cands):>9.0f}"
    )

allus = [num(r["us_exact"]) for r in rows]
allus = [u for u in allus if u is not None]
print(f"{'ALL':<10} {len(allus):>5} {sum(allus) / 1e6:>9.3f} {st.median(allus):>10.1f} "
      f"{pct(allus, 90):>10.1f} {pct(allus, 99):>10.1f} {max(allus):>11.1f}")

# speedup of the pipeline against the exact stage alone
tot_native = sum(num(r["us_native"]) or 0.0 for r in rows)
tot_exact = sum(allus)
print(f"\npipeline total {tot_native / 1e6:.3f}s   exact stage alone {tot_exact / 1e6:.3f}s "
      f"= {tot_exact / tot_native:.2f}x the pipeline total over the whole corpus")

# accuracy of the exact stage against the answer the pipeline returns
short, same, over = [], 0, []
for r in rows:
    re_, rn = num(r["radius_exact"]), num(r["radius_native"])
    if re_ is None or rn is None or rn == 0:
        continue
    rel = (rn - re_) / abs(rn)
    if abs(rel) < 1e-9:
        same += 1
    elif rel > 0:
        short.append(rel)
    else:
        over.append(rel)
print("\n=== exact stage against the returned answer ===")
print(f"identical to 1e-9 relative: {same}")
if short:
    print(f"short of the returned radius: {len(short)} cases, median {st.median(short):.3e}, "
          f"p90 {pct(short, 90):.3e}, max {max(short):.3e}")
    for thr in (1e-6, 1e-4, 1e-2, 1e-1):
        print(f"  short by more than {thr:.0e}: {sum(1 for s in short if s > thr)}")
if over:
    print(f"larger than the returned radius: {len(over)} cases, max {max(over):.3e}")

# the exact circle must fit, measured against the polygon
bad = []
for r in rows:
    re_, ce = num(r["radius_exact"]), num(r["clearance_exact"])
    if re_ is None or ce is None or ce == 0:
        continue
    rel = (re_ - ce) / ce
    if rel > 0:
        bad.append((rel, r["name"]))
bad.sort(reverse=True)
print(f"exact radii above their own measured clearance: {len(bad)}")
for cls, thr in (("rounding level (<=1e-12)", 1e-12), ("small (<=1e-9)", 1e-9), ("real (>1e-9)", 1e-9)):
    got = [b for b in bad if b[0] <= thr] if thr == 1e-12 else [b for b in bad if b[0] <= thr and b[0] > 1e-12] if thr == 1e-9 else [b for b in bad if b[0] > 1e-9]
    print(f"  {cls}: {len(got)}")
if bad:
    print("  worst relative overclaims:")
    for rel, name in bad[:5]:
        print(f"    {rel:.3e}  {name}")

# scaling on real data
print("\n=== exact stage by vertex count (real polygons only) ===")
buckets = [(0, 10), (10, 25), (25, 50), (50, 100), (100, 250), (250, 1000), (1000, 10**9)]
real = by_kind.get("real", [])
print(f"{'verts':<14} {'n':>5} {'median us':>10} {'max us':>11} {'median cand':>12}")
for lo, hi in buckets:
    rs = [r for r in real if lo <= (num(r["verts"]) or 0) < hi]
    us = [num(r["us_exact"]) for r in rs]
    us = [u for u in us if u is not None]
    if not us:
        continue
    cands = [num(r["candidates_exact"]) or 0 for r in rs]
    print(f"{f'{lo}-{hi}':<14} {len(us):>5} {st.median(us):>10.1f} {max(us):>11.1f} {st.median(cands):>12.0f}")
