"""02_tilt_scan.py - is the linear end-on geometry a minimum? (env ml)

01_bind_n2.py found that a 25-deg tilted start relaxes lower in energy than
the linear starts (which can stay linear by symmetry: no sideways force).
Here: starts at Li-N 2.1 A tilted 15 and 30 deg in 8 azimuths around the
clearance axis, N2 only relaxed (same calculator). Lowest -> complex_relaxed.extxyz.
"""
import json
from pathlib import Path

import numpy as np
from ase.io import read, write
from ase.optimize import BFGS

HERE = Path(__file__).parent
from ase import Atoms
from ase.constraints import FixAtoms
import torch
from mace.calculators import mace_mp

ROOT = HERE.parents[1]
prev = json.loads((HERE / "binding.json").read_text())
LI = prev["Li_index"]
fw = read(ROOT / "structure" / "LiLSX_P1_MACE.cif")
A = fw.cell[0, 0]
CENTRE = np.full(3, 3 / 8) * A
li = fw.positions[LI]
lin = read(HERE / "starts" / "d2.2_relaxed.extxyz")
AXIS = lin.positions[-1] - lin.positions[-2]
AXIS /= np.linalg.norm(AXIS)
e1 = np.cross(AXIS, [0, 0, 1.0]); e1 /= np.linalg.norm(e1)
e2 = np.cross(AXIS, e1)

torch.cuda.set_per_process_memory_fraction(4.0 / 7.96)
calc = mace_mp(model="medium", device="cuda", default_dtype="float32", dispersion=True)


def geom(at):
    n1, n2 = len(at) - 2, len(at) - 1
    if at.get_distance(LI, n2, mic=True) < at.get_distance(LI, n1, mic=True):
        n1, n2 = n2, n1
    return {"Li-N_A": round(float(at.get_distance(LI, n1, mic=True)), 3),
            "Li-N-N_deg": round(float(at.get_angle(LI, n1, n2, mic=True)), 2),
            "N-N_A": round(float(at.get_distance(n1, n2, mic=True)), 4)}


rows = []
for tilt in (15, 30):
    for az in range(0, 360, 45):
        t, p = np.radians(tilt), np.radians(az)
        d = np.cos(t) * AXIS + np.sin(t) * (np.cos(p) * e1 + np.sin(p) * e2)
        n1 = li + 2.1 * AXIS
        at = fw.copy() + Atoms("N2", positions=[n1, n1 + 1.098 * d])
        at.set_constraint(FixAtoms(indices=list(range(len(fw)))))
        at.calc = calc
        BFGS(at, logfile=None, maxstep=0.1).run(fmax=0.005, steps=600)
        g = geom(at)
        g.update({"tilt0": tilt, "azimuth0": az, "E_eV": float(at.get_potential_energy())})
        rows.append(g)
        write(HERE / "starts" / f"tilt{tilt}_az{az}_relaxed.extxyz", at)
        print(g)
emin = min(r["E_eV"] for r in rows + [{"E_eV": v["E_eV"]} for v in prev["runs"].values()])
for r in rows:
    r["dE_meV"] = round(1000 * (r["E_eV"] - emin), 1)
best = min(rows, key=lambda r: r["E_eV"])
lin_E = prev["runs"]["d2.2"]["E_eV"]
summary = {"scan": rows, "best": best,
           "linear_minus_best_meV": round(1000 * (lin_E - best["E_eV"]), 1),
           "linear_minus_best_kJmol": round(96.485 * (lin_E - best["E_eV"]), 2)}
(HERE / "tilt_scan.json").write_text(json.dumps(summary, indent=1))
src = HERE / "starts" / f"tilt{best['tilt0']}_az{best['azimuth0']}_relaxed.extxyz"
if best["E_eV"] < min(v["E_eV"] for v in prev["runs"].values()):
    write(HERE / "complex_relaxed.extxyz", read(src))
print(json.dumps({k: v for k, v in summary.items() if k != "scan"}, indent=1))
