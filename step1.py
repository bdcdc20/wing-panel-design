"""
ME32012 Wing Panel Design - Part (i), Step 1: element geometry
SAAB 340, Group C. Elements numbered 1 (root) to 7 (tip), tabulated tip -> root.
"""
import pandas as pd

# ---- Aircraft data (brief p.2) ----
S = 40.58            # total wing area incl. fuselage rectangle (m^2)
y_body = 0.8         # body side from centreline (m)
c_body = 2.7         # chord at body side (m)
c_tip = 1.0          # tip chord (m)

# Element lengths from Table 1, element 1 (root) -> 7 (tip)
lengths = [1.3, 1.0, 1.3, 1.3, 1.3, 1.8, 1.8]

y_tip = y_body + sum(lengths)            # 10.6 m
taper = (c_body - c_tip) / (y_tip - y_body)


def chord(y):
    """Local chord (m) at distance y (m) from the centreline."""
    if y <= y_body:
        return c_body                     # constant through fuselage
    return c_body - taper * (y - y_body)


# ---- Build geometry table ----
rows = []
y_in = y_body
for i, b in enumerate(lengths, start=1):
    y_out = y_in + b
    c_in, c_out = chord(y_in), chord(y_out)
    c_mean = (c_in + c_out) / 2
    rows.append(dict(element=i, y_in=y_in, y_out=y_out, b=b,
                     c_in=c_in, c_out=c_out, c_mean=c_mean,
                     A=b * c_mean))
    y_in = y_out

geo = pd.DataFrame(rows).set_index("element").sort_index(ascending=False)  # tip -> root

# ---- Checks ----
A_wing = geo["A"].sum()
A_fus = y_body * c_body
print(geo.round(3).to_string())
print(f"\nSemi-span (tip from centreline): {y_tip:.2f} m")
print(f"Sum of element areas (one wing): {A_wing:.3f} m^2   (expect 18.13)")
print(f"+ fuselage rectangle {A_fus:.2f} m^2 = {A_wing + A_fus:.3f} m^2   (expect S/2 = {S/2:.2f})")
assert abs(A_wing + A_fus - S / 2) < 1e-6, "Geometry does not reproduce S/2"
