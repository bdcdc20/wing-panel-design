"""
ME32012 Wing Panel Design - Part (i), Step 2: limit SF, BM and torque (reference axes)
SAAB 340, Group E (element 6). Max-g manoeuvre at Vc, MTOW, strip method (Butler pp. 11-15).
Values are accumulated tip -> root and quoted at each element's INNER boundary.
Sign convention: shear/bending positive for net upward load; torque positive leading edge up.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from step1 import geo, S, chord, y_tip

# ---- Data (brief p.2) ----
W = 130.37                 # MTOW (kN)
W_lb = W * 1e3 / 4.44822   # MTOW (lb) for CS 25.337
n = max(2.1 + 24000 / (W_lb + 10000), 2.5)   # limit load factor, = 2.71
Vc = 140.0                 # m/s EAS
rho0 = 1.225               # sea-level density (EAS => rho0)
q = 0.5 * rho0 * Vc**2 / 1e3   # dynamic pressure (kN/m^2)
CM0 = -0.08
g = 9.81
rho_fuel = 900.0           # kg/m^3 (SG 0.9)
W_eng = 4.5                # engine + main LG on element 2 (kN)
x_eng = 1.0                # engine c.g. ahead of reference axis (m)
fuel_total = 25.8          # kN, both wings
tanks = [1, 3, 4, 5]       # fuelled elements (2, 6, 7 empty)

# Reference axis (40% c) = wing-box centre (20-60% c). Lift at 25% c is 0.15c ahead.
# Wing weight and fuel act on the reference axis -> no pitching moment.

# ---- Fuel: tank capacity = box between spars, 0.4c x 0.12c x b (mean chord) ----
geo = geo.copy()
geo["tank_cap"] = 0.0
geo.loc[tanks, "tank_cap"] = (rho_fuel * g / 1e3) * (0.4 * geo["c_mean"]) * (0.12 * geo["c_mean"]) * geo["b"]
cap_wing = geo["tank_cap"].sum()
fuel_wing = fuel_total / 2
# Fill outboard tanks first (inboard empties first). Capacity (~10.1 kN) < 12.9 kN per wing,
# so all four tanks are full; the excess is taken as carried outside the wing (fuselage),
# which is conservative for wing bending (less inertia relief).
remaining = fuel_wing
geo["fuel"] = 0.0
for e in sorted(tanks, reverse=True):
    geo.loc[e, "fuel"] = min(geo.loc[e, "tank_cap"], remaining)
    remaining -= geo.loc[e, "fuel"]

# ---- Element loads at limit (n g) ----
geo["L"] = n * (W / S) * geo["A"]
geo["Ww"] = n * 0.15 * W * geo["A"] / S
geo["Wf"] = n * geo["fuel"]
geo["We"] = 0.0
geo.loc[2, "We"] = n * W_eng
geo["Fs"] = geo["L"] - geo["Ww"] - geo["Wf"] - geo["We"]
geo["Mps"] = (geo["L"] * 0.15 * geo["c_mean"]
              + CM0 * q * geo["b"] * geo["c_mean"]**2
              - geo["We"] * x_eng)

# ---- Accumulate tip -> root (Omega = 0: M_REF = M_R, T_REF = M_PS) ----
SF, BM, T = 0.0, 0.0, 0.0
rows = [dict(element="tip", y=y_tip, SF=0.0, BM=0.0, T=0.0)]
for e, r in geo.iterrows():           # geo is ordered tip -> root
    BM += SF * r["b"] + r["Fs"] * r["b"] / 2
    SF += r["Fs"]
    T += r["Mps"]
    rows.append(dict(element=e, y=r["y_in"], SF=SF, BM=BM, T=T))
loads = pd.DataFrame(rows)

# ---- Design station: inboard edge of element 6 ----
y6 = geo.loc[6, "y_in"]
c6 = chord(y6)
M6 = loads.loc[loads["element"] == 6, "BM"].item()
# Check: exact integration of the uniform-pressure + wing-weight load (no point loads outboard of 6)
k = n * (W / S) * 0.85
yy = np.linspace(y6, y_tip, 20001)
M6_exact = np.trapezoid(k * np.array([chord(y) for y in yy]) * (yy - y6), yy)
Nx = lambda M: 1.5 * M * 1e6 / ((0.12 * c6 * 1e3) * (0.4 * c6 * 1e3))   # N/mm, ultimate

if __name__ == "__main__":
    pd.set_option("display.width", 140)
    print(f"n = {n:.3f}   q = {q:.3f} kN/m^2")
    print(f"Tank capacity per wing = {cap_wing:.2f} kN vs fuel per wing {fuel_wing:.2f} kN "
          f"-> {remaining:.2f} kN carried outside wing")
    print("\nElement loads (limit, kN / kNm):")
    print(geo[["y_in", "y_out", "A", "L", "Ww", "Wf", "We", "Fs", "Mps"]].round(3).to_string())
    print("\nLimit SF / BM / T at element inner boundaries (reference axes):")
    print(loads.round(3).to_string(index=False))
    print(f"\nDesign station y = {y6:.1f} m, c = {c6:.3f} m, box {0.4*c6*1e3:.0f} x {0.12*c6*1e3:.0f} mm")
    print(f"  Strip method  M = {M6:.2f} kNm  -> ultimate Nx = {Nx(M6):.1f} N/mm")
    print(f"  Exact integr. M = {M6_exact:.2f} kNm  -> ultimate Nx = {Nx(M6_exact):.1f} N/mm")

    geo.to_csv("loads_elements.csv", float_format="%.4f")
    loads.to_csv("loads_sf_bm_t.csv", index=False, float_format="%.4f")

    fig, axes = plt.subplots(3, 1, figsize=(7, 9), sharex=True)
    for ax, col, lab in zip(axes, ["SF", "BM", "T"],
                            ["Shear force (kN)", "Bending moment (kNm)", "Torque (kNm)"]):
        ax.plot(loads["y"], loads[col], "o-", color="k", ms=4)
        for yb in loads["y"]:
            ax.axvline(yb, color="0.85", lw=0.6, zorder=0)
        ax.axvspan(y6, geo.loc[6, "y_out"], color="tab:blue", alpha=0.08)
        ax.axhline(0, color="0.5", lw=0.6)
        ax.set_ylabel(lab)
        ax.grid(False)
    axes[0].set_title(f"Limit loads, n = {n:.2f}, MTOW, Vc = {Vc:.0f} m/s EAS (element 6 shaded)")
    axes[-1].set_xlabel("Spanwise position from aircraft centreline, y (m)")
    axes[-1].set_xlim(0.8, y_tip)
    fig.tight_layout()
    fig.savefig("fig_sf_bm_t.png", dpi=200)
    print("\nSaved loads_elements.csv, loads_sf_bm_t.csv, fig_sf_bm_t.png")
