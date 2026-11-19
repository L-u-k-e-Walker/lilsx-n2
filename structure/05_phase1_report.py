"""05_phase1_report.py - charge balance, cell, atom counts, Li sites and
supercage contents of LiLSX_P1_MACE.cif. Writes results/phase1_report.json.

Formal charges: Si +4, Al +3, O -2, Li +1 (framework [AlSiO4]- per Al).
Supercage centres are the 8b-type points (3/8,3/8,3/8) of Fd-3m origin 2.
"""
import csv
import json
from pathlib import Path

import numpy as np
from ase.io import read
from ase.spacegroup import crystal

HERE = Path(__file__).parent
at = read(HERE / "LiLSX_P1_MACE.cif")
sym = np.array(at.get_chemical_symbols())
counts = {e: int((sym == e).sum()) for e in ("Si", "Al", "O", "Li")}
q = {"Si": 4, "Al": 3, "O": -2, "Li": 1}
framework_charge = sum(q[e] * counts[e] for e in ("Si", "Al", "O"))
total_charge = framework_charge + counts["Li"]

sites = {}
with open(HERE / "results" / "li_sites_P1.csv") as fh:
    for r in csv.DictReader(fh):
        sites[int(r["atom_index_0based"])] = r["site"]
assert all(sym[i] == "Li" for i in sites) and len(sites) == counts["Li"]

# supercage centres; the relaxed cell is cubic, so fractional = conventional
sc = crystal("X", [(3/8, 3/8, 3/8)], spacegroup=227, setting=2,
             cellpar=[1, 1, 1, 90, 90, 90]).get_scaled_positions()
assert len(sc) == 8
f = at.get_scaled_positions()
cell = at.cell[:]
per_cage = []
for k, c in enumerate(sc):
    d = f - c
    d -= np.round(d)
    r = np.linalg.norm(d @ cell, axis=1)
    near = [i for i in sites if r[i] < 8.0]
    per_cage.append({"cage": k, "centre_frac": np.round(c, 4).tolist(),
                     "Li_within_8A": {s: sum(sites[i] == s for i in near) for s in ("SI'", "SII", "SIII")},
                     "Li_distance_from_centre_A": {s: sorted(round(float(r[i]), 2) for i in near if sites[i] == s)
                                                   for s in ("SII", "SIII")}})

rep = {
    "file": "structure/LiLSX_P1_MACE.cif",
    "cell_A": np.round(at.cell.lengths(), 4).tolist(),
    "cell_angles_deg": np.round(at.cell.angles(), 4).tolist(),
    "volume_A3": round(float(at.get_volume()), 2),
    "atoms_total": len(at), "counts": counts,
    "Si/Al": counts["Si"] / counts["Al"],
    "charge_framework_e": framework_charge, "charge_Li_e": counts["Li"],
    "charge_total_e": total_charge,
    "Li_sites": {s: sum(v == s for v in sites.values()) for s in ("SI'", "SII", "SIII")},
    "supercages": per_cage,
}
(HERE / "results" / "phase1_report.json").write_text(json.dumps(rep, indent=1))
print(json.dumps({k: v for k, v in rep.items() if k != "supercages"}, indent=1))
for c in per_cage:
    print(c)
assert total_charge == 0
