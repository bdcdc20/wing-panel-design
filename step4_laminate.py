"""
ME32012 Wing Panel Design - Part (ii), Step 4: discrete laminate design and stacking sequences
SAAB 340, Group E (element 6). Turns the continuous optimum of step 3 into buildable laminates
of 0.125 mm plies, re-optimises b_st for the discrete thicknesses and re-checks every constraint.

Ply-count rules (Butler p.25):
  - angles 0, +45, -45, 90 only; balanced and symmetric -> +/-45 plies in multiples of 4;
  - 0 and 90 counts even, except that ONE of them may be odd (a single ply on the mid-plane);
  - at least 10 % of each angle (0, +45, -45, 90);
  - percentages within PCT_TOL of the target breakdown (44/44/12 skin, 60/30/10 stiffener);
    7 points is the nearest achievable with +/-45 in multiples of 4 at these thin gauges.
Moduli are recomputed with Krenchel from the ACTUAL ply counts.
Stacking-sequence rules checked on the full laminate: symmetric, balanced, <= 3 adjacent
plies of the same angle, +/-45 pair at both surfaces, every +45 next to a -45, >= 10 % each.
"""
import itertools
import json
import math
import numpy as np

import step3_panel as s3
from safe_save import save_json
from step3_panel import props, feasible, E11, bsk, m_Al, Nx

T_PLY = 0.125                     # mm
PCT_TOL = 0.07
TARGET = dict(skin=np.array([0.44, 0.44, 0.12]), stiffener=np.array([0.60, 0.30, 0.10]))
BST_GRID = np.arange(30.0, 60.01, 0.5)


def laminate(n0, n45, n90):
    N = n0 + n45 + n90
    return dict(n0=n0, n45=n45, n90=n90, N=N, t=N * T_PLY,
                pct=np.array([n0, n45, n90]) / N,
                E=E11 * (n0 + 0.25 * n45) / N)          # Krenchel: cos^4 45 = 0.25, cos^4 90 = 0


def candidates(target, t_max, tol=PCT_TOL, odd_mid=True):
    out = []
    for n0 in range(1, int(t_max / T_PLY) + 1):
        for n45 in range(4, int(t_max / T_PLY) + 1, 4):
            for n90 in range(1, int(t_max / T_PLY) + 1):
                if (n0 % 2) + (n90 % 2) > (1 if odd_mid else 0):
                    continue
                lam = laminate(n0, n45, n90)
                p = lam["pct"]
                if lam["t"] > t_max or p[0] < 0.1 or p[1] < 0.2 or p[2] < 0.1:
                    continue                       # p[1] >= 0.2 -> +45 and -45 each >= 10 %
                if np.max(np.abs(p - target)) > tol:
                    continue
                out.append(lam)
    return out


def best_bst(sk, st):
    """Lightest feasible b_st on the grid for a given skin/stiffener laminate pair."""
    p = props(BST_GRID, sk["t"], st["t"], Esk=sk["E"], Est=st["E"])
    m = np.where(feasible(p), p["mass"], np.inf)
    i = int(np.argmin(m))
    return (m[i], BST_GRID[i]) if np.isfinite(m[i]) else (np.inf, None)


def search(sk_list, st_list):
    res = []
    for sk, st in itertools.product(sk_list, st_list):
        m, b = best_bst(sk, st)
        if np.isfinite(m):
            res.append((m, b, sk, st))
    return sorted(res, key=lambda r: r[0])


# ---------------- stacking sequences ----------------
def check_stack(stack):
    """Return list of rule violations for a full ply list (strings '0', '+45', '-45', '90')."""
    bad = []
    N = len(stack)
    if any(a not in ("0", "+45", "-45", "90") for a in stack):
        bad.append("angle")
    if stack != stack[::-1]:
        bad.append("symmetric")
    if stack.count("+45") != stack.count("-45"):
        bad.append("balanced")
    run = max(len(list(g)) for _, g in itertools.groupby(stack))
    if run > 3:
        bad.append(f"grouping {run}")
    if set(stack[:2]) != {"+45", "-45"} or set(stack[-2:]) != {"+45", "-45"}:
        bad.append("surface +/-45")
    for i, a in enumerate(stack):
        if a in ("+45", "-45"):
            opp = "-45" if a == "+45" else "+45"
            nb = stack[max(i - 1, 0):i] + stack[i + 1:i + 2]
            if opp not in nb:
                bad.append("+/-45 adjacent")
                break
    for a in ("0", "+45", "-45", "90"):
        if stack.count(a) < 0.1 * N:
            bad.append(f">=10% {a}")
    return bad


def stacking_sequence(lam):
    """Depth-first search for a half-stack that, mirrored, satisfies every rule.
    Blocks: P = (+45, -45) pair, '0', '90'. Choices ordered to spread plies evenly."""
    h = {"P": lam["n45"] // 4, "0": lam["n0"] // 2, "90": lam["n90"] // 2}
    mid = "0" if lam["n0"] % 2 else ("90" if lam["n90"] % 2 else None)
    total = sum(h.values())

    def expand(blocks):
        plies = []
        for b in blocks:
            plies += ["+45", "-45"] if b == "P" else [b]
        return plies

    def full(blocks):
        half = expand(blocks)
        return half + ([mid] if mid else []) + half[::-1]

    def dfs(blocks, left):
        if not any(left.values()):
            st = full(blocks)
            return st if not check_stack(st) else None
        plies = expand(blocks)
        # prefer the block type that lags its share the most (spreads 0s, 90s and pairs)
        order = sorted((k for k in left if left[k] > 0),
                       key=lambda k: (h[k] - left[k]) / h[k] - len(blocks) / total)
        for k in order:
            if k != "P":
                run = 0
                for a in reversed(plies):
                    if a != k:
                        break
                    run += 1
                if run >= 3:
                    continue
            left[k] -= 1
            r = dfs(blocks + [k], left)
            left[k] += 1
            if r:
                return r
        return None

    left = dict(h)
    left["P"] -= 1                               # start with a +/-45 pair at the surface
    return dfs(["P"], left)


SUB = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")


def notation(stack):
    """Compact notation of the upper half, e.g. [±45/0₂/90/...]s, mid-plane ply overbarred."""
    n = len(stack)
    half, mid = stack[: n // 2], (stack[n // 2] if n % 2 else None)
    parts, i = [], 0
    while i < len(half):
        if half[i] == "+45" and i + 1 < len(half) and half[i + 1] == "-45":
            k = 0
            while half[i:i + 2] == ["+45", "-45"]:
                k += 1
                i += 2
            parts.append("±45" + (str(k).translate(SUB) if k > 1 else ""))
            continue
        k = 1
        while i + k < len(half) and half[i + k] == half[i]:
            k += 1
        parts.append(half[i] + (str(k).translate(SUB) if k > 1 else ""))
        i += k
    if mid:
        parts.append("".join(ch + "̅" for ch in mid))  # overbar = single mid-plane ply
    return "[" + "/".join(parts) + "]s"


# ---------------- run ----------------
opt = s3.opt                                     # continuous optimum from step 3

# (a) baseline: round each orientation of the step-3 thicknesses up (even counts, +/-45 x4)
def round_up(t, pct):
    n = [t * p / T_PLY for p in pct]
    return laminate(2 * math.ceil(n[0] / 2), 4 * math.ceil(n[1] / 4), 2 * math.ceil(n[2] / 2))

base_sk = round_up(opt["tsk"], TARGET["skin"])
base_st = round_up(opt["tst"], TARGET["stiffener"])
base_m, base_b = best_bst(base_sk, base_st)

# (b) search of rule-compliant discrete laminates near the target breakdowns
sk_list = candidates(TARGET["skin"], 4.0)
st_list = candidates(TARGET["stiffener"], 2.6)
ranked = search(sk_list, st_list)
m_fin, b_fin, sk, st = ranked[0]

# (c) sensitivity: looser percentage tolerance
loose = search(candidates(TARGET["skin"], 4.0, tol=0.10), candidates(TARGET["stiffener"], 2.6, tol=0.10))

p = props(b_fin, sk["t"], st["t"], Esk=sk["E"], Est=st["E"])
seq_sk, seq_st = stacking_sequence(sk), stacking_sequence(st)
assert seq_sk and not check_stack(seq_sk) and seq_st and not check_stack(seq_st)
assert feasible(p)


def lam_str(l):
    return (f"{l['N']} plies ({l['n0']}/{l['n45']}/{l['n90']}), t = {l['t']:.3f} mm, "
            f"{'/'.join(f'{100*x:.0f}' for x in l['pct'])} %, E = {l['E']/1e3:.2f} kN/mm²")


if __name__ == "__main__":
    print(f"Nx = {Nx:.1f} N/mm; continuous optimum (step 3): b_st {opt['bst']:.1f}, "
          f"t_sk {opt['tsk']:.3f}, t_st {opt['tst']:.3f} mm, {opt['mass']:.3f} kg/m²\n")

    print("(a) Round each orientation up (even counts, ±45 in 4s):")
    print("    skin     ", lam_str(base_sk))
    print("    stiffener", lam_str(base_st))
    print(f"    -> best b_st {base_b:.1f} mm, {base_m:.3f} kg/m²\n")

    print(f"(b) Search: {len(sk_list)} skin x {len(st_list)} stiffener laminates "
          f"(±{100*PCT_TOL:.0f} pts of target, one mid-plane ply allowed). Lightest five:")
    for m, b, a, c in ranked[:5]:
        print(f"    {m:6.3f} kg/m²  b_st {b:4.1f}  skin {a['n0']}/{a['n45']}/{a['n90']} "
              f"({a['t']:.3f})  stiff {c['n0']}/{c['n45']}/{c['n90']} ({c['t']:.3f})")
    m, b, a, c = loose[0]
    print(f"\n(c) Sensitivity, ±10 pts: {m:.3f} kg/m² at b_st {b:.1f}, skin "
          f"{a['n0']}/{a['n45']}/{a['n90']} {'/'.join(f'{100*x:.0f}' for x in a['pct'])} %, "
          f"stiff {c['n0']}/{c['n45']}/{c['n90']} {'/'.join(f'{100*x:.0f}' for x in c['pct'])} %")

    print("\n================ FINAL DESIGN ================")
    print(f"b_st = {b_fin:.1f} mm, b_sk = {bsk:.0f} mm, flanges 2 x 30 mm, rib pitch 600 mm")
    print("Skin     ", lam_str(sk))
    print("          ", notation(seq_sk), "  plies:", " ".join(seq_sk))
    print("Stiffener", lam_str(st))
    print("          ", notation(seq_st), "  plies:", " ".join(seq_st))
    print(f"A_st = {p['Ast']:.1f} mm², Et = {p['Et']/1e3:.1f} kN/mm, z̄ = {p['zbar']:.2f} mm, "
          f"EI = {p['EI']/1e9:.2f}e6 kN·mm²")
    print(f"\n{'Mode':28s}{'Critical Nx (N/mm)':>20s}{'RF':>8s}{'Required':>10s}")
    for k, name, req in [("op", "Skin between stiffeners", 1.0), ("cl", "Skin under stiffener", 1.0),
                         ("st", "Stiffener local", 1.0), ("mat", "Material (strain 0.0045)", 1.0),
                         ("eu", "Euler (overall)", 1.5)]:
        rf = p[f"RF_{k}"]
        print(f"{name:28s}{rf*Nx:20.0f}{rf:8.2f}{req:10.2f}")
    print(f"Applied ultimate Nx = {Nx:.0f} N/mm")
    print(f"\nMass = {p['mass']:.3f} kg/m²  (continuous {opt['mass']:.3f}; rounding penalty "
          f"{100*(p['mass']/opt['mass']-1):.1f} %)")
    print(f"Saving vs aluminium {m_Al} kg/m²: {100*(1-p['mass']/m_Al):.1f} %")

    save_json(dict(bst=b_fin, bsk=bsk, bf=s3.bf, tsk=sk["t"], tst=st["t"],
                   skin=dict(n0=sk["n0"], n45=sk["n45"], n90=sk["n90"], E=sk["E"], stack=seq_sk,
                             notation=notation(seq_sk)),
                   stiffener=dict(n0=st["n0"], n45=st["n45"], n90=st["n90"], E=st["E"],
                                  stack=seq_st, notation=notation(seq_st)),
                   mass=p["mass"], RF={k: p[f"RF_{k}"] for k in ["op", "cl", "st", "mat", "eu"]},
                   Nx=Nx), "final_design.json", indent=2, ensure_ascii=False)
    print("\nSaved final_design.json")
