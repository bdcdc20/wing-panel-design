"""
ME32012 Wing Panel Design - Part (ii), Step 3: minimum-mass CFRP hat-stiffened panel
SAAB 340, Group E: compression panel AC at the inboard edge of element 6 (y = 7.0 m).
Continuous optimum (before ply rounding). Method: Butler notes pp. 21-23 + tips sheet:
  1. choose b_st;
  2. equate skin and stiffener local buckling strain, using the larger of the
     open (b_sk - 2b_st - b_f) and closed (2b_st) skin breadths;
  3. scale thicknesses until the most critical local mode carries N_x exactly,
     then check material strain (RF >= 1) and Euler buckling (RF >= 1.5);
  4. repeat over b_st and take the minimum of m = rho (t_sk + A_st / b_sk).
Units: N, mm (N/mm^2 = MPa).
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from safe_save import save_csv, save_fig
from step2_loads import M6, Nx as Nx_from_M, c6

# ---- Data ----
E11 = 130e3                                   # N/mm^2
Esk = E11 * (0.44 + 0.44 * 0.25 + 0.12 * 0)   # Krenchel, 44/44/12 -> 71.5 kN/mm^2
Est = E11 * (0.60 + 0.30 * 0.25 + 0.10 * 0)   # Krenchel, 60/30/10 -> 87.75 kN/mm^2
Ksk, Kst = 2.23, 1.65
eps_max = 0.0045
rho = 1600e-9                                 # kg/mm^3
bsk = 200.0                                   # stiffener pitch (mm)
bf = 30.0                                     # width of each flange (mm)
L_rib = 600.0                                 # rib spacing, element 6 (mm)
RF_euler_min = 1.5
m_Al = 9.1                                    # aluminium panel, element 6 (kg/m^2)

Nx = Nx_from_M(M6)                            # ultimate load per unit width (N/mm), ~701
r3 = np.sqrt(3.0)


def props(bst, tsk, tst, bsk=bsk, bf=bf, L=L_rib, Nx=Nx, Esk=Esk, Est=Est):
    """Stiffness, mass and reserve factors of one repeating portion (vectorised)."""
    Ast = tst * (3 * bst + 2 * bf)            # hat top + two webs (length b_st) + two flanges
    Et = Esk * tsk + Est * Ast / bsk          # smeared axial stiffness (N/mm)

    # Neutral axis from skin mid-plane (Butler p.22); webs at 60 deg -> hat height (sqrt3/2) b_st
    zf = (tsk + tst) / 2                      # flanges
    zw = zf + r3 / 4 * bst                    # web centroids
    zt = zf + r3 / 2 * bst                    # hat top
    EA_sk, EA_f = Esk * bsk * tsk, Est * 2 * bf * tst
    EA_w, EA_t = Est * 2 * bst * tst, Est * bst * tst
    zb = (EA_f * zf + EA_w * zw + EA_t * zt) / (EA_sk + EA_f + EA_w + EA_t)

    # Bending stiffness about the neutral axis (Butler p.23); each web EI_cg = Est t b^3 / 16
    EI = (EA_sk * zb**2 + Esk * bsk * tsk**3 / 12
          + EA_f * (zf - zb)**2 + Est * 2 * bf * tst**3 / 12
          + EA_t * (zt - zb)**2 + Est * bst * tst**3 / 12
          + 2 * Est * tst * bst**3 / 16 + EA_w * (zw - zb)**2)

    b_op = bsk - 2 * bst - bf                 # between flange centres
    b_cl = 2 * bst                            # under the hat
    N = dict(op=Et * Ksk * (tsk / b_op)**2,
             cl=Et * Ksk * (tsk / b_cl)**2,
             st=Et * Kst * (tst / bst)**2,
             mat=Et * eps_max,
             eu=np.pi**2 * EI / (L**2 * bsk))
    out = {f"RF_{k}": v / Nx for k, v in N.items()}
    out.update(Et=Et, Ast=Ast, EI=EI, zbar=zb,
               mass=rho * (tsk + Ast / bsk) * 1e6)          # kg/m^2
    return out


def feasible(p):
    return (np.minimum.reduce([p["RF_op"], p["RF_cl"], p["RF_st"], p["RF_mat"]]) >= 1.0) \
        & (p["RF_eu"] >= RF_euler_min)


def design(bst):
    """Minimum-thickness design for a given b_st with equal skin/stiffener buckling strain."""
    b_max = max(bsk - 2 * bst - bf, 2 * bst)
    ratio = (bst / b_max) * np.sqrt(Ksk / Kst)   # t_st / t_sk from Ksk(tsk/bmax)^2 = Kst(tst/bst)^2
    lo, hi = 0.01, 50.0
    for _ in range(100):                         # bisection on t_sk (all RFs increase with t)
        t = 0.5 * (lo + hi)
        if feasible(props(bst, t, ratio * t)):
            hi = t
        else:
            lo = t
    tsk = hi
    p = props(bst, tsk, ratio * tsk)
    limits = dict(op=1.0, cl=1.0, st=1.0, mat=1.0, eu=RF_euler_min)
    active = "+".join(k for k, lim in limits.items() if p[f"RF_{k}"] < lim * (1 + 1e-3))
    return dict(bst=bst, tsk=tsk, tst=ratio * tsk, active=active, **p)


# ---- 1. Validate against Butler's worked example (notes pp. 21-23) ----
ex = props(45.0, 3.0, 2.0, bsk=180.0, bf=35.0, L=800.0, Nx=1500.0, Esk=72e3, Est=88e3)
ex_RF = [ex[k] for k in ["RF_op", "RF_cl", "RF_st", "RF_mat", "RF_eu"]]
# Local and material RFs match to 2 d.p.; z_bar = 9.44 mm matches exactly. Butler quotes
# EI = 16.4e6 but his own p.23 expression evaluates to 16.7e6 here, so Euler RF = 0.95 (vs 0.94).
assert abs(ex["zbar"] - 9.44) < 0.01, ex["zbar"]
assert np.allclose(ex_RF[:4], [1.84, 0.69, 0.91, 1.25], atol=0.01), ex_RF
assert abs(ex_RF[4] - 0.94) < 0.02, ex_RF

# ---- 2. Sweep b_st (open breadth must stay positive: b_st < (b_sk - b_f)/2 = 85 mm) ----
sweep = pd.DataFrame([design(b) for b in np.arange(15.0, 80.01, 0.1)])
opt = sweep.loc[sweep["mass"].idxmin()]

# ---- 3. Independent check: brute-force grid over (t_sk, t_st) without the equal-strain rule ----
TSK, TST = np.meshgrid(np.arange(1.0, 5.0, 0.005), np.arange(0.5, 4.0, 0.005), indexing="ij")
best_grid = (np.inf, None)
for b in np.arange(30.0, 60.01, 0.5):
    p = props(b, TSK, TST)
    m = np.where(feasible(p), p["mass"], np.inf)
    i = np.unravel_index(np.argmin(m), m.shape)
    if m[i] < best_grid[0]:
        best_grid = (m[i], (b, TSK[i], TST[i]))

if __name__ == "__main__":
    pd.set_option("display.width", 140)
    print(f"Design station y = 7.0 m, c = {c6:.3f} m; M_limit = {M6:.2f} kNm")
    print(f"Ultimate Nx = 1.5 M / (0.12c x 0.4c) = {Nx:.1f} N/mm")
    print(f"Esk = {Esk/1e3:.2f} kN/mm^2, Est = {Est/1e3:.2f} kN/mm^2")
    print("Butler worked-example RFs reproduced:", np.round(ex_RF, 2))

    print("\nSweep (every 5 mm):")
    cols = ["bst", "tsk", "tst", "mass", "RF_op", "RF_cl", "RF_st", "RF_mat", "RF_eu", "active"]
    print(sweep[np.isclose(sweep["bst"] % 5, 0) | np.isclose(sweep["bst"] % 5, 5)][cols]
          .round(3).to_string(index=False))

    print("\nContinuous optimum:")
    for k in cols:
        v = opt[k]
        print(f"  {k:7s} = {v:.3f}" if isinstance(v, float) else f"  {k:7s} = {v}")
    print(f"  A_st    = {opt['Ast']:.1f} mm^2, smeared t = {opt['tsk'] + opt['Ast']/bsk:.3f} mm")
    print(f"  Mass saving vs aluminium {m_Al} kg/m^2: {100*(1-opt['mass']/m_Al):.1f} %")
    mg, (bg, tg, sg) = best_grid
    print(f"\nBrute-force check (free t_sk, t_st): m = {mg:.3f} kg/m^2 at "
          f"b_st = {bg:.1f}, t_sk = {tg:.3f}, t_st = {sg:.3f}")

    save_csv(sweep, "panel_sweep.csv", index=False, float_format="%.4f")

    # ---- Figure: mass and reserve factors vs b_st (two panels, shared x, one y-scale each) ----
    blue, orange, ink, muted = "#2a78d6", "#eb6834", "#222222", "#8a8a85"
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7, 7), sharex=True,
                                   gridspec_kw=dict(height_ratios=[1.1, 1]))
    x = sweep["bst"]
    eu_active = sweep["active"].str.contains("eu")
    for ax in (ax1, ax2):
        ax.axvline(opt["bst"], color=muted, lw=1, ls="--", zorder=1)
        if eu_active.any():
            ax.axvspan(x[eu_active].min(), x[eu_active].max(), color="0.93", zorder=0)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)

    ax1.plot(x, sweep["mass"], color=ink, lw=2)
    ax1.plot(opt["bst"], opt["mass"], "o", color=ink, ms=8, mec="white", mew=2, zorder=5)
    ax1.annotate(f"optimum  b$_{{st}}$ = {opt['bst']:.1f} mm, {opt['mass']:.2f} kg/m²",
                 (opt["bst"], opt["mass"]), xytext=(8, -18), textcoords="offset points",
                 color=ink, fontsize=9)
    ax1.axhline(m_Al, color=muted, lw=1, ls=":")
    ax1.annotate("aluminium panel 9.1 kg/m²", (44, m_Al), xytext=(0, 4),
                 textcoords="offset points", ha="left", va="bottom", color=muted, fontsize=9)
    ax1.set_ylabel("Mass per unit area (kg/m²)")
    ax1.set_ylim(bottom=0)
    if eu_active.any():
        ax1.text(x[eu_active].mean(), 0.5, "Euler\nactive", ha="center", va="bottom",
                 color=muted, fontsize=8)

    ax2.plot(x, sweep["RF_eu"], color=blue, lw=2, label="Euler buckling RF")
    ax2.plot(x, sweep["RF_mat"], color=orange, lw=2, label="Damage-tolerant strain RF")
    ax2.axhline(RF_euler_min, color=blue, lw=1, ls=":")
    ax2.axhline(1.0, color=orange, lw=1, ls=":")
    ax2.annotate("Euler RF ≥ 1.5", (x.max(), RF_euler_min), xytext=(0, 3), textcoords="offset points",
                 ha="right", va="bottom", fontsize=8, color=muted)
    ax2.annotate("strain RF ≥ 1.0", (x.max(), 1.0), xytext=(0, -3), textcoords="offset points",
                 ha="right", va="top", fontsize=8, color=muted)
    ax2.set_ylabel("Reserve factor")
    ax2.set_xlabel("Stiffener width b$_{st}$ (mm)")
    ax2.set_ylim(0, 8)                        # Euler RF keeps rising beyond (17 at b_st = 80 mm)
    ax2.legend(frameon=False, loc="upper right", bbox_to_anchor=(1.0, 0.70), fontsize=9)

    ax1.set_title(f"Element 6 hat-stiffened CFRP panel: N$_x$ = {Nx:.0f} N/mm (ultimate),\n"
                  f"b$_{{sk}}$ = 200 mm, rib pitch 600 mm, equal skin/stiffener buckling strain",
                  fontsize=10)
    fig.tight_layout()
    save_fig(fig, "fig_panel_optimisation.png", dpi=200)
    print("\nSaved panel_sweep.csv, fig_panel_optimisation.png")
