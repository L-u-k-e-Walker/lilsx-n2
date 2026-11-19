"""01_bind_n2.py - N2 end-on at an exposed SIII Li+ in Li-LSX (env ml).

1. Choose the SIII Li in the supercage centred at (3/8,3/8,3/8): prefer 2-fold
   O-coordinated Li (most exposed), seen side-on by the Phase-4 camera, which
   looks along +[111] with the front half of the cage cut away.
2. Find the N2 axis direction from that Li with the largest clearance from the
   framework (cone search around the Li -> cage-centre vector).
3. Relax only the two N atoms (framework and all Li fixed) with MACE-MP-0
   medium + D3(BJ), float32 (same settings as structure/02_relax.py), from
   several starts: Li-N 2.2 / 2.6 / 3.0 A on the clearance axis, and 2.4 A
   tilted 25 deg. MACE is used only to get a plausible geometry.
4. Single points for a rough MACE interaction energy (not a validated number).

Outputs: complex_relaxed.extxyz (lowest-energy start), binding.json, starts/*.extxyz
"""
import json
from pathlib import Path

import numpy as np
import torch
from ase import Atoms
from ase.constraints import FixAtoms
from ase.io import read, write
from ase.optimize import BFGS
from mace.calculators import mace_mp

HERE = Path(__file__).parent
ROOT = HERE.parents[1]
(HERE / "starts").mkdir(exist_ok=True)

fw = read(ROOT / "structure" / "LiLSX_P1_MACE.cif")
A = fw.cell[0, 0]
CENTRE = np.full(3, 3 / 8) * A
VIEW = np.ones(3) / np.sqrt(3)            # Phase-4 camera direction (camera -> centre)
site = {}
for line in (ROOT / "structure" / "results" / "li_sites_P1.csv").read_text().splitlines()[1:]:
    p = line.split(",")
    site[int(p[0])] = (p[1], int(p[5]))   # (site label, n O within 2.5 A)
sym = np.array(fw.get_chemical_symbols())


def mic(v):
    return v - np.round(v / A) * A


# ---------------------------------------------------------------- choose Li
cands = []
for i, (lab, cn) in site.items():
    if lab != "SIII":
        continue
    v = mic(fw.positions[i] - CENTRE)
    r = np.linalg.norm(v)
    if r > 8.0:
        continue
    u = v / r
    depth = float(u @ VIEW)               # >0: behind the centre (kept by the cutaway)
    cands.append({"index": i, "cn": cn, "r_centre": round(float(r), 3), "view_depth": round(depth, 3)})
# exposed (CN 2) first, then side-on to the camera: depth near 0.3
cands.sort(key=lambda c: (c["cn"], abs(c["view_depth"] - 0.3)))
LI = cands[0]["index"]
li = fw.positions[LI]
print("SIII Li in this supercage:", cands)
print("chosen Li index", LI)

# ---------------------------------------------------------------- clearance axis
others = np.array([j for j in range(len(fw)) if j != LI])


def clearance(axis, d1=2.4, nn=1.098):
    pts = np.array([li + d1 * axis, li + (d1 + nn) * axis])
    dd = mic(fw.positions[others][None, :, :] - pts[:, None, :])
    return float(np.linalg.norm(dd, axis=2).min())


to_c = mic(CENTRE - li)
to_c /= np.linalg.norm(to_c)
rng = np.random.default_rng(1)
best = (clearance(to_c), to_c)
for _ in range(4000):
    v = rng.normal(size=3)
    v /= np.linalg.norm(v)
    ax = to_c + 1.2 * v
    ax /= np.linalg.norm(ax)
    if np.degrees(np.arccos(np.clip(ax @ to_c, -1, 1))) > 60:
        continue
    c = clearance(ax)
    if c > best[0]:
        best = (c, ax)
AXIS = best[1]
print(f"clearance axis: {np.degrees(np.arccos(AXIS @ to_c)):.1f} deg from Li->centre, "
      f"min N-framework distance at Li-N 2.4 A: {best[0]:.2f} A")


def perp(a):
    t = np.cross(a, [0, 0, 1.0])
    if np.linalg.norm(t) < 1e-3:
        t = np.cross(a, [0, 1.0, 0])
    return t / np.linalg.norm(t)


def make(d_lin, tilt_deg=0.0, nn=1.098):
    t = np.radians(tilt_deg)
    n_dir = np.cos(t) * AXIS + np.sin(t) * perp(AXIS)
    n1 = li + d_lin * AXIS
    n2 = n1 + nn * n_dir
    at = fw.copy() + Atoms("N2", positions=[n1, n2])
    at.set_constraint(FixAtoms(indices=list(range(len(fw)))))
    return at


def geom(at):
    n1, n2 = len(at) - 2, len(at) - 1
    d1 = at.get_distance(LI, n1, mic=True)
    d2 = at.get_distance(LI, n2, mic=True)
    if d2 < d1:
        n1, n2, d1 = n2, n1, d2
    ang = at.get_angle(LI, n1, n2, mic=True)
    nn = at.get_distance(n1, n2, mic=True)
    dn = at.get_all_distances(mic=True)[[n1, n2]][:, :len(fw)]
    dn[0, LI] = dn[1, LI] = 99
    near = int(np.argmin(dn.min(axis=0)))
    return {"Li-N_A": round(float(d1), 3), "Li-N-N_deg": round(float(ang), 2), "N-N_A": round(float(nn), 4),
            "closest_other_atom": f"{sym[near]}{near}", "closest_other_A": round(float(dn.min()), 3),
            "N2_axis_vs_Li_to_centre_deg": round(float(np.degrees(np.arccos(
                np.clip((at.positions[n2] - at.positions[n1]) / nn @ to_c, -1, 1)))), 1)}


# ---------------------------------------------------------------- relax N2 only
torch.cuda.set_per_process_memory_fraction(4.0 / 7.96)
calc = mace_mp(model="medium", device="cuda", default_dtype="float32", dispersion=True)
starts = {"d2.2": (2.2, 0), "d2.6": (2.6, 0), "d3.0": (3.0, 0), "d2.4_tilt25": (2.4, 25)}
runs = {}
for name, (d, tilt) in starts.items():
    at = make(d, tilt)
    at.calc = calc
    g0 = geom(at)
    opt = BFGS(at, logfile=str(HERE / "starts" / f"{name}.log"), maxstep=0.1)
    conv = opt.run(fmax=0.005, steps=600)
    runs[name] = {"start": g0, "final": geom(at), "E_eV": float(at.get_potential_energy()),
                  "converged": bool(conv), "steps": opt.nsteps}
    write(HERE / "starts" / f"{name}_relaxed.extxyz", at)
    print(name, runs[name])

best_name = min(runs, key=lambda k: runs[k]["E_eV"])
final = read(HERE / "starts" / f"{best_name}_relaxed.extxyz")
write(HERE / "complex_relaxed.extxyz", final)

# ---------------------------------------------------------------- rough energy
fw_only = fw.copy()
fw_only.calc = calc
e_fw = fw_only.get_potential_energy()
n2 = Atoms("N2", positions=[[0, 0, 0], [0, 0, 1.098]], cell=[20, 20, 20], pbc=True)
n2.calc = calc
BFGS(n2, logfile=None).run(fmax=0.001)
e_n2 = n2.get_potential_energy()
d_gas = n2.get_distance(0, 1)
# N2 far from any Li (cage centre) for a same-cell reference
far = fw.copy() + Atoms("N2", positions=[CENTRE - [0, 0, 0.549], CENTRE + [0, 0, 0.549]])
far.calc = calc
e_far = far.get_potential_energy()

out = {
    "Li_index": LI, "Li_site": "SIII", "Li_O_coordination": site[LI][1],
    "Li_candidates_in_cage": cands,
    "axis_vs_Li_to_centre_deg": round(float(np.degrees(np.arccos(AXIS @ to_c))), 1),
    "runs": runs, "best_start": best_name, "geometry": runs[best_name]["final"],
    "N2_gas_bond_MACE_A": round(float(d_gas), 4),
    "E_int_vs_gas_N2_kJmol": round((runs[best_name]["E_eV"] - e_fw - e_n2) * 96.485, 1),
    "E_int_vs_N2_at_cage_centre_kJmol": round((runs[best_name]["E_eV"] - e_far) * 96.485, 1),
    "note": "MACE-MP-0 + D3 numbers, framework rigid; for a plausible geometry only",
}
(HERE / "binding.json").write_text(json.dumps(out, indent=1))
print(json.dumps({k: v for k, v in out.items() if k not in ("runs", "Li_candidates_in_cage")}, indent=1))
