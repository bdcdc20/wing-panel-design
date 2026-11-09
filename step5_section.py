"""
ME32012 Wing Panel Design - Part (ii), Step 5: to-scale cross-section of the final design
SAAB 340, Group E (element 6). Reads final_design.json (step 4) and draws
  (a) one repeating portion of the hat-stiffened panel, true scale, dimensioned (cf. brief Fig. 2);
  (b) the skin and stiffener ply stacks (stacking sequences, ply thickness 0.125 mm).
Geometry: skin on top; hat bonded to its underside by two 30 mm flanges; webs of length b_st
at 60 deg; hat top width b_st; hat centreline height (sqrt3/2) b_st (as in Butler p.22).
"""
import json
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Rectangle

from safe_save import save_fig

d = json.load(open("final_design.json"))
bsk, bst, bf = d["bsk"], d["bst"], d["bf"]
tsk, tst = d["tsk"], d["tst"]
h = np.sqrt(3) / 2 * bst                      # hat centreline height
dx = bst / 2                                  # web horizontal projection (b_st cos 60)

INK, MUTED, FILL = "#222222", "#6f6f6a", "#cfd3d8"
ANGLE_STYLE = {                               # colour + hatch (secondary encoding) + label
    "0":   dict(color="#2a78d6", hatch="---"),
    "+45": dict(color="#eb6834", hatch="///"),
    "-45": dict(color="#1baf7a", hatch="\\\\\\"),
    "90":  dict(color="#eda100", hatch="|||"),
}


def offset_polyline(pts, t):
    """Closed polygon of a polyline thickened by t (centred), mitred joints."""
    pts = np.asarray(pts, float)
    seg = np.diff(pts, axis=0)
    nrm = np.c_[-seg[:, 1], seg[:, 0]] / np.linalg.norm(seg, axis=1)[:, None]
    side = []
    for s in (+1, -1):
        out = [pts[0] + s * t / 2 * nrm[0]]
        for i in range(1, len(pts) - 1):
            m = nrm[i - 1] + nrm[i]
            m /= np.linalg.norm(m)
            out.append(pts[i] + s * t / 2 * m / np.dot(m, nrm[i]))
        out.append(pts[-1] + s * t / 2 * nrm[-1])
        side.append(np.array(out))
    return np.vstack([side[0], side[1][::-1]])


def dim(ax, p1, p2, text, off=(0, 0), txt_off=(0, 0), ha="center", va="center", ext=True):
    """Dimension line between p1 and p2 shifted by off, with extension lines."""
    p1, p2, off = np.array(p1, float), np.array(p2, float), np.array(off, float)
    a, b = p1 + off, p2 + off
    ax.annotate("", a, b, arrowprops=dict(arrowstyle="<|-|>", color=MUTED, lw=0.8,
                                          shrinkA=0, shrinkB=0, mutation_scale=7))
    if ext and np.any(off):
        for p, q in ((p1, a), (p2, b)):
            ax.plot(*zip(p, q + 0.15 * off / max(np.linalg.norm(off), 1e-9) * 3),
                    color=MUTED, lw=0.5)
    m = (a + b) / 2 + np.array(txt_off)
    ax.text(*m, text, ha=ha, va=va, fontsize=8, color=INK,
            bbox=dict(fc="white", ec="none", pad=0.6))


# ---- Hat centreline (z = 0 at skin/flange interface, z negative downwards) ----
xc = bsk / 2
zf = -tst / 2                                 # flange mid-plane
x_root = xc - bst / 2 - dx                    # web roots: closed width = 2 b_st
hat = [(x_root - bf, zf), (x_root, zf), (x_root + dx, zf - h),
       (xc + bst / 2, zf - h), (xc + bst / 2 + dx, zf), (xc + bst / 2 + dx + bf, zf)]

fig = plt.figure(figsize=(11, 5.2))
gs = fig.add_gridspec(1, 2, width_ratios=[3.3, 1], wspace=0.08)
ax = fig.add_subplot(gs[0])
ax2 = fig.add_subplot(gs[1])

# ---- (a) cross-section ----
ax.add_patch(Rectangle((0, 0), bsk, tsk, fc=FILL, ec=INK, lw=0.8))
ax.add_patch(Polygon(offset_polyline(hat, tst), closed=True, fc=FILL, ec=INK, lw=0.8))
ax.plot([0, 0], [-tst - 9, tsk + 6], color=MUTED, lw=0.6, ls=(0, (6, 3)))
ax.plot([bsk, bsk], [-tst - 9, tsk + 6], color=MUTED, lw=0.6, ls=(0, (6, 3)))

dim(ax, (0, tsk), (bsk, tsk), f"b$_{{sk}}$ = {bsk:.0f}", off=(0, 20))
dim(ax, (x_root, tsk), (2 * xc - x_root, tsk), f"closed skin 2b$_{{st}}$ = {2*bst:.0f}",
    off=(0, 8))
fl_c = x_root - bf / 2
dim(ax, (x_root - bf, -tst), (x_root, -tst), f"{bf:.0f}", off=(0, -5))
dim(ax, (xc + bst / 2 + dx, -tst), (xc + bst / 2 + dx + bf, -tst), f"{bf:.0f}", off=(0, -5))
dim(ax, (xc - bst / 2, zf - h - tst / 2), (xc + bst / 2, zf - h - tst / 2),
    f"b$_{{st}}$ = {bst:.1f}", off=(0, -6))
# web length, parallel to the right-hand web
w0, w1 = np.array(hat[4]), np.array(hat[3])
n = np.array([np.sin(np.radians(60)), np.cos(np.radians(60))])  # outward normal of right web
dim(ax, w1, w0, f"b$_{{st}}$ = {bst:.1f}", off=-7 * n, txt_off=(-2, 0), ha="right", ext=False)
# hat centreline height, inside the hat (left of centre, clear of the web dimension)
dim(ax, (xc - 14, zf), (xc - 14, zf - h), f"{h:.1f}", txt_off=(-1.5, 0), ha="right", ext=False)
# 60 deg between flange/skin and web (inside the hat at the left web root)
arc = np.radians(np.linspace(0, 60, 20))
ax.plot(x_root + 11 * np.cos(arc), zf - 11 * np.sin(arc), color=MUTED, lw=0.7)
ax.text(x_root + 13, zf - 3.5, "60°", fontsize=8, color=INK, va="top")
# neighbouring repeat (ghost): skin and right-hand flange, for the open-skin width
ax.add_patch(Rectangle((-70, 0), 70, tsk, fc="none", ec=MUTED, lw=0.6, ls=(0, (3, 2))))
ax.add_patch(Rectangle((x_root + bst + 2 * dx - bsk, -tst), bf, tst, fc="none", ec=MUTED,
                       lw=0.6, ls=(0, (3, 2))))
nb_c = x_root + bst + 2 * dx - bsk + bf / 2      # neighbour flange centre
ax.plot([nb_c, nb_c], [-tst, -tst - 18], color=MUTED, lw=0.5)
ax.plot([fl_c, fl_c], [-tst, -tst - 18], color=MUTED, lw=0.5)
dim(ax, (nb_c, -tst - 16), (fl_c, -tst - 16), "", ext=False)
ax.text((nb_c + fl_c) / 2, -tst - 19, f"open skin b$_{{sk}}$−2b$_{{st}}$−30 = {bsk-2*bst-bf:.0f}\n"
        "(between flange centres)", fontsize=7.5, color=INK, ha="center", va="top")
ax.text(-35, tsk + 3, "adjacent repeat", fontsize=7, color=MUTED, ha="center", va="bottom")

# thickness call-outs
ax.annotate(f"skin t$_{{sk}}$ = {tsk:.3f} mm\n{d['skin']['n0']+d['skin']['n45']+d['skin']['n90']} plies "
            f"{d['skin']['notation']}", (172, tsk), xytext=(150, 34), fontsize=8, color=INK,
            arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.6))
ax.annotate(f"hat stiffener + flanges t$_{{st}}$ = {tst:.2f} mm\n"
            f"{d['stiffener']['n0']+d['stiffener']['n45']+d['stiffener']['n90']} plies "
            f"{d['stiffener']['notation']}", (xc + bst / 2 - 4, zf - h - tst / 2),
            xytext=(138, -50), va="top", fontsize=8, color=INK,
            arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.6))

ax.set_xlim(-72, bsk + 40)
ax.set_ylim(-h - 34, tsk + 50)
ax.set_aspect("equal")
ax.axis("off")
ax.set_title("(a) Repeating portion, true scale (mm)", fontsize=10, loc="left")

# ---- (b) ply stacks ----
def draw_stack(ax, x0, stack, title):
    n = len(stack)
    for i, a in enumerate(stack):                # ply 1 = top surface
        y = n - 1 - i
        st = ANGLE_STYLE[a]
        ax.add_patch(Rectangle((x0, y), 1.0, 1.0, fc=st["color"], ec="white", lw=1.0,
                               hatch=st["hatch"]))
        ax.text(x0 + 1.12, y + 0.5, a + "°", fontsize=7, va="center", color=INK)
    # mid-plane (through the middle of the single centre ply when n is odd)
    ax.plot([x0 - 0.15, x0 + 1.15], [n / 2, n / 2], color=INK, lw=0.6, ls=":")
    ax.text(x0 + 0.5, n + 0.6, title, ha="center", va="bottom", fontsize=8, color=INK)

plt.rcParams["hatch.color"] = (1, 1, 1, 0.55)
plt.rcParams["hatch.linewidth"] = 0.6
draw_stack(ax2, 0.0, d["skin"]["stack"], f"Skin\n{len(d['skin']['stack'])} plies, {tsk:.3f} mm")
draw_stack(ax2, 2.6, d["stiffener"]["stack"],
           f"Stiffener\n{len(d['stiffener']['stack'])} plies, {tst:.2f} mm")
ax2.text(1.8, -0.9, "dotted line = mid-plane (symmetric)", ha="center", fontsize=7, color=MUTED)
ax2.set_xlim(-0.3, 4.3)
ax2.set_ylim(-1.5, len(d["skin"]["stack"]) + 3)
ax2.axis("off")
ax2.set_title("(b) Stacking sequences (0.125 mm plies,\nnot to scale through thickness)",
              fontsize=10, loc="left")

fig.suptitle(f"Element 6 CFRP hat-stiffened panel: final design, {d['mass']:.2f} kg/m² "
             f"(aluminium 9.1 kg/m²)", fontsize=11, x=0.02, ha="left")
save_fig(fig, "fig_cross_section.png", dpi=300, bbox_inches="tight")
save_fig(fig, "fig_cross_section.svg", bbox_inches="tight")
print("Saved fig_cross_section.png, fig_cross_section.svg")
