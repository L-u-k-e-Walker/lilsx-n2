"""02_relax.py - relax Li-LSX candidates with MACE-MP-0 (medium) + D3(BJ).

usage: python 02_relax.py <candidate id> [fmax]

Stage 1: atoms only, fixed starting cell (a = 24.714 A), fmax 0.05 eV/A.
Stage 2: atoms + cubic cell (hydrostatic strain), fmax (default 0.01 eV/A)
         -> build/LiLSX_c<id>_relaxed.extxyz (main result)
Stage 3: atoms + free cell from stage 2 -> build/LiLSX_c<id>_relaxed_freecell.extxyz
Per-candidate summary: build/relax_c<id>.json.
float32 model (memory, see below). Run in env `ml` (mace-torch 0.3.17, ASE 3.29), same setup as hkust1_mace/.
"""
import json
import sys
import time
from pathlib import Path

import numpy as np
import torch
from ase.filters import FrechetCellFilter
from ase.io import read, write
from ase.optimize import FIRE, BFGS
from mace.calculators import mace_mp

HERE = Path(__file__).parent
B = HERE / "build"
cid = int(sys.argv[1])
fmax = float(sys.argv[2]) if len(sys.argv) > 2 else 0.01

atoms = read(B / f"LiLSX_init_c{cid}.extxyz")
# float32: peak 3.3 GB reserved for 672 atoms vs 6.6 GB in float64 (8 GB card,
# WSL spills VRAM to RAM instead of failing); energy differs by 2e-4 eV.
torch.cuda.set_per_process_memory_fraction(4.0 / 7.96)
atoms.calc = mace_mp(model="medium", device="cuda", default_dtype="float32",
                     dispersion=True)
t0 = time.time()
e0 = atoms.get_potential_energy()

# stage 1: FIRE is robust far from the minimum (Li start positions are rough)
FIRE(atoms, logfile=str(B / f"relax_c{cid}_stage1.log"),
     trajectory=str(B / f"relax_c{cid}_stage1.traj")).run(fmax=0.05, steps=1500)
e1 = atoms.get_potential_energy()

# Cell DOF scaling: FrechetCellFilter converts stress to "forces" as
# virial / exp_cell_factor (default exp_cell_factor = N atoms = 672, which lets
# ~0.1 GPa through at fmax 0.02). exp_cell_factor = 50 with fmax 0.01 eV/A means
# |stress| < ~0.005 GPa at convergence.
CELL_FACTOR = 50.0


def relax_cell(atoms, tag, **kw):
    opt = BFGS(FrechetCellFilter(atoms, exp_cell_factor=CELL_FACTOR, **kw),
               logfile=str(B / f"relax_c{cid}_{tag}.log"),
               trajectory=str(B / f"relax_c{cid}_{tag}.traj"), maxstep=0.1)
    return opt.run(fmax=fmax, steps=3000)


def summary_of(atoms, conv):
    return {"converged": bool(conv), "E_eV": atoms.get_potential_energy(),
            "cell_lengths_A": atoms.cell.lengths().tolist(),
            "cell_angles_deg": atoms.cell.angles().tolist(),
            "volume_A3": atoms.get_volume(),
            "max_force_eV_A": float(np.linalg.norm(atoms.get_forces(), axis=1).max()),
            "stress_GPa": (atoms.get_stress() * 160.21766).tolist()}


# stage 2: cubic cell (hydrostatic strain only) + atoms. Main result: the
# experimental room-temperature Li-LSX cell is cubic (SIII disorder averages out)
conv2 = relax_cell(atoms, "stage2_cubic", hydrostatic_strain=True)
cubic = summary_of(atoms, conv2)
write(B / f"LiLSX_c{cid}_relaxed.extxyz", atoms)

# stage 3: free cell from the cubic minimum, to measure how far the ordered
# SIII arrangement wants to distort the cell
conv3 = relax_cell(atoms, "stage3_free")
free = summary_of(atoms, conv3)
write(B / f"LiLSX_c{cid}_relaxed_freecell.extxyz", atoms)

summary = {"candidate": cid, "fmax": fmax, "cell_factor": CELL_FACTOR,
           "E_start_eV": e0, "E_fixedcell_eV": e1, "converged": cubic["converged"],
           "E_relaxed_eV": cubic["E_eV"], "cell_lengths_A": cubic["cell_lengths_A"],
           "cell_angles_deg": cubic["cell_angles_deg"], "volume_A3": cubic["volume_A3"],
           "max_force_eV_A": cubic["max_force_eV_A"], "cubic": cubic, "free_cell": free,
           "wall_s": time.time() - t0}
(B / f"relax_c{cid}.json").write_text(json.dumps(summary, indent=1))
print(json.dumps(summary, indent=1))
