import pymupdf
import numpy as np
from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import haloconcentration as hc


PDF_PATH = "tmp/arxiv_2408.06515v2_source/Prada.pdf"
MASSES = np.array([1.0e3, 1.0e6, 1.0e9, 1.0e12])

page = pymupdf.open(PDF_PATH)[0]
curve_drawings = [
    drawing
    for drawing in page.get_drawings()
    if drawing.get("width") == 1.5 and len(drawing.get("items", [])) == 49
]

# Prada.pdf 的坐标轴由 Matplotlib 直接保存为矢量。横轴是 z=0--12；纵轴
# 通过两个相邻主刻度 C=5 和 C=10 建立线性坐标变换。
x_left = 59.32993698120117
x_right = 433.1357116699219
y_at_c5 = 262.8352966308594
y_at_c10 = 217.80650329589844

source_curves = []
for drawing in curve_drawings:
    items = drawing["items"]
    points = [items[0][1]] + [item[2] for item in items]
    x_pdf = np.array([point.x for point in points])
    y_pdf = np.array([point.y for point in points])
    redshift = 12.0 * (x_pdf - x_left) / (x_right - x_left)
    concentration = 5.0 + 5.0 * (y_at_c5 - y_pdf) / (y_at_c5 - y_at_c10)
    source_curves.append(concentration)

source_curves = np.asarray(source_curves)

old_curves = []
new_curves = []
for mass0 in MASSES:
    old_mass_track = hc.mass_accretion_history(
        mass0,
        redshift,
        "prada12",
    )
    new_mass_track = hc.mass_accretion_history(
        mass0,
        redshift,
        "prada12_hmf_sigma",
    )
    old_curves.append(
        hc.concentration_prada12(
            old_mass_track,
            redshift,
            cap_high_peak=True,
        )
    )
    new_curves.append(
        hc.concentration_prada12_hmf_sigma(
            new_mass_track,
            redshift,
            cap_high_peak=True,
        )
    )

old_curves = np.asarray(old_curves)
new_curves = np.asarray(new_curves)

for label, curves in [("old Eq.23", old_curves), ("new HMF sigma", new_curves)]:
    relative = (curves - source_curves) / source_curves
    print(label)
    print("  MAPE all points =", np.mean(np.abs(relative)))
    print("  RMS relative    =", np.sqrt(np.mean(relative**2)))
    print("  max relative    =", np.max(np.abs(relative)))
    print("  z=0 values      =", curves[:, 0])

print("source z=0 values =", source_curves[:, 0])
print("redshift endpoints =", redshift[0], redshift[-1], "point count", redshift.size)
print("SOURCE_REDSHIFT =", repr(redshift.tolist()))
print("SOURCE_PRADA_CURVES =", repr(source_curves.tolist()))
