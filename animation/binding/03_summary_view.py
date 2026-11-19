"""03_summary_view.py - final bound geometry numbers + a check render (env viz).

Reads complex_relaxed.extxyz (lowest-energy N2-only relaxation), binding.json
and tilt_scan.json; writes binding_summary.json and bound_view.png (same
cutaway camera as structure/06_view_supercage.py, active Li highlighted).
"""
import warnings
warnings.filterwarnings("ignore", message=".*OVITO.*PyPI")
import json
from pathlib import Path

import numpy as np
from ase.io import read, write

HERE = Path(__file__).parent
b = json.loads((HERE / "binding.json").read_text())
t = json.loads((HERE / "tilt_scan.json").read_text())
at = read(HERE / "complex_relaxed.extxyz")
LI = b["Li_index"]
nfw = len(at) - 2
A = at.cell[0, 0]
CENTRE = np.full(3, 3 / 8) * A
n1, n2 = nfw, nfw + 1
if at.get_distance(LI, n2, mic=True) < at.get_distance(LI, n1, mic=True):
    n1, n2 = n2, n1
D = at.get_all_distances(mic=True)
sym = np.array(at.get_chemical_symbols())
fwmask = np.arange(len(at)) < nfw
fwmask[LI] = False


def mic(v):
    return v - np.round(v / A) * A


ax = mic(at.positions[n2] - at.positions[n1])
ax /= np.linalg.norm(ax)
to_c = mic(CENTRE - at.positions[LI])
to_c /= np.linalg.norm(to_c)
near = {}
for n, lab in ((n1, "N_inner"), (n2, "N_outer")):
    j = np.where(fwmask)[0][np.argsort(D[n][fwmask])[:3]]
    near[lab] = [f"{sym[k]}{k}: {D[n, k]:.2f} A" for k in j]
lio = np.sort(D[LI][sym == "O"])[:3]
E_best = t["best"]["E_eV"]
e_gas_ref = b["runs"][b["best_start"]]["E_eV"] - b["E_int_vs_gas_N2_kJmol"] / 96.485  # E_fw + E_N2
summary = {
    "Li_index": LI, "Li_site": "SIII",
    "Li-N_A": round(float(D[LI, n1]), 3),
    "Li-N-N_deg": round(float(at.get_angle(LI, n1, n2, mic=True)), 2),
    "N-N_A": round(float(D[n1, n2]), 4), "N-N_gas_MACE_A": b["N2_gas_bond_MACE_A"],
    "Li-N2_centre_of_mass_A": round(float(np.linalg.norm(mic(0.5 * (at.positions[n1] + at.positions[n2]) - at.positions[LI]))), 3),
    "N2_axis_vs_Li_to_cage_centre_deg": round(float(np.degrees(np.arccos(ax @ to_c))), 1),
    "Li-O_framework_3_nearest_A": np.round(lio, 3).tolist(),
    "nearest_framework_contacts": near,
    "E_int_vs_gas_N2_kJmol_MACE": round((E_best - e_gas_ref) * 96.485, 1),
    "linear_end_on_above_minimum_kJmol": t["linear_minus_best_kJmol"],
    "second_basin_169deg_above_minimum_kJmol": round(33.8e-3 * 96.485, 1),
}
(HERE / "binding_summary.json").write_text(json.dumps(summary, indent=1))
print(json.dumps(summary, indent=1))

# ---------------------------------------------------------------- render
from ovito.io.ase import ase_to_ovito
from ovito.pipeline import Pipeline, StaticSource
from ovito.modifiers import CreateBondsModifier, DeleteSelectedModifier
from ovito.vis import Viewport, TachyonRenderer, BondsVis

VIEW = np.ones(3) / np.sqrt(3)
d = mic(at.positions - CENTRE)
keep = (np.linalg.norm(d, axis=1) <= 9.5) & ((d @ VIEW) >= -1.5)
keep[[n1, n2]] = True
sub = at[np.where(keep)[0]]
sub.positions = CENTRE + d[keep]
idx = np.where(keep)[0]
data = ase_to_ovito(sub)
data.cell_.vis.enabled = False
col = data.particles_.create_property("Color")
rad = data.particles_.create_property("Radius")
for k, i in enumerate(idx):
    s = sym[i]
    col[k], rad[k] = {"O": ((0.80, 0.35, 0.30), 0.18), "Si": ((0.75, 0.70, 0.55), 0.22),
                      "Al": ((0.60, 0.70, 0.80), 0.22), "Li": ((0.55, 0.25, 0.80), 0.70),
                      "N": ((0.15, 0.45, 1.00), 0.55)}[s]
    if i == LI:
        col[k], rad[k] = (0.80, 0.45, 1.00), 0.85
p = Pipeline(source=StaticSource(data=data))
cb = CreateBondsModifier(mode=CreateBondsModifier.Mode.Pairwise)
for tt in ("Si", "Al"):
    cb.set_pairwise_cutoff(tt, "O", 2.0)
cb.set_pairwise_cutoff("N", "N", 1.3)
cb.vis.width = 0.22
cb.vis.coloring_mode = BondsVis.ColoringMode.ByParticle
p.modifiers.append(cb)
p.add_to_scene()
vp = Viewport(type=Viewport.Type.Perspective, camera_pos=tuple(CENTRE - 45 * VIEW), camera_dir=tuple(VIEW),
              camera_up=(0, 0, 1), fov=np.radians(28))
vp.render_image(size=(960, 720), filename=str(HERE / "bound_view.png"), background=(1, 1, 1),
                renderer=TachyonRenderer(ambient_occlusion=True, antialiasing_samples=12))
print("wrote bound_view.png")
