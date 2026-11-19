"""01_build.py - build an unrelaxed Li-LSX (Li96 Al96 Si96 O384) cell.

Inputs (sources/):
  FAU_IZA.cif  IZA-SC FAU framework (pure SiO2, DLS-optimised, a = 24.345 A,
               Fd-3m origin choice 2)
  1530083.cif  COD, Herden et al., Zeolites 2 (1982) 131, Li-X (Si/Al 1.23),
               used for the experimental SI' and SII Li x-coordinates
  1529943.cif  COD, Forano et al., J. Solid State Chem. 82 (1989) 95, anhydrous
               Li-Y neutron structure; second reference for SI'/SII and the cell

Steps
  1. Expand the IZA framework to P1 (192 T, 384 O).
  2. Si/Al = 1 ordering: the FAU T-T graph is bipartite, so a two-colouring
     gives strict Si-O-Al alternation (Loewenstein). Checked below.
  3. Rescale the cell to a_start (Li-X, Herden 1982) so the framework starts
     near the right size; MACE relaxes the cell later.
  4. Li sites (96 per cell, one per Al):
       SI'  32e (x,x,x), sodalite side of the D6R 6-ring, x from Herden (LiX)
       SII  32e (x,x,x), supercage side of the sodalite 6-ring, x from Herden
       SIII 48f (x,1/8,1/8), on the -4 axis in front of the sodalite 4-ring,
            facing the supercage; 32 of the 48 sites are filled (occ. 2/3)
     Occupancies SI' 32 / SII 32 / SIII 32 follow the dehydrated Li-LSX
     literature (Plevert et al. 1997; Feuerstein & Lobo 1998): SI' and SII are
     full, and the remaining 32 Li sit in the supercage at SIII/III'.
  5. Choose which 32 of the 48 SIII sites are filled by Monte Carlo on a
     Li-Li Coulomb sum (Ewald). All 48f sites are symmetry-equivalent in Fd-3,
     so only the Li-Li arrangement differs between choices. The NCAND lowest
     distinct arrangements are written for MACE screening in 02_relax.py.
"""
import json
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.io import read, write
from ase.spacegroup import crystal

HERE = Path(__file__).parent
OUT = HERE / "build"
OUT.mkdir(exist_ok=True)
rng = np.random.default_rng(20261118)

A_START = 24.714          # A, Li-X (Si/Al 1.23), Herden 1982 (COD 1530083)
X_SI1P = 0.0392           # SI'  x,x,x  (Herden 1982)
X_SII = 0.2275            # SII  x,x,x  (Herden 1982)
D_SIII = 2.0              # A, start: Li this far beyond the 4-ring centre
NCAND = 3

# ---------------------------------------------------------------- framework
fau = read(HERE / "sources" / "FAU_IZA.cif")
assert len(fau) == 576
fau.set_cell(np.eye(3) * A_START, scale_atoms=True)
cell = fau.cell[:]
frac = fau.get_scaled_positions()
isT = fau.numbers == 14
Tidx = np.where(isT)[0]
Oidx = np.where(~isT)[0]


def mic(d):
    d = d - np.round(d)
    return d @ cell


# T-O bonds (< 2.0 A), each O must bridge exactly two T
nbr = {}
for o in Oidx:
    r = np.linalg.norm(mic(frac[Tidx] - frac[o]), axis=1)
    t = Tidx[r < 2.0]
    assert len(t) == 2, (o, len(t))
    nbr[o] = t
adj = {t: [] for t in Tidx}
for t1, t2 in nbr.values():
    adj[t1].append(t2)
    adj[t2].append(t1)

# bipartite two-colouring by BFS
colour = {Tidx[0]: 0}
queue = [Tidx[0]]
while queue:
    t = queue.pop()
    for u in adj[t]:
        if u not in colour:
            colour[u] = 1 - colour[t]
            queue.append(u)
        else:
            assert colour[u] != colour[t], "T graph not bipartite"
assert len(colour) == 192
symbols = np.array(fau.get_chemical_symbols(), dtype=object)
for t, c in colour.items():
    symbols[t] = "Al" if c else "Si"
n_al = int((symbols == "Al").sum())
assert n_al == 96
# Loewenstein: every O bridges one Si and one Al
assert all({symbols[a], symbols[b]} == {"Si", "Al"} for a, b in nbr.values())

# ---------------------------------------------------------------- Li sites
# Fd-3m origin choice 2 (same setting as the IZA and COD files):
#   16c (0,0,0) = D6R centre, 8a (1/8,1/8,1/8) = sodalite centre,
#   8b (3/8,3/8,3/8) = supercage centre (checked in the asserts below).
#   32e (x,x,x): SI' for x < 1/8 near the D6R, SII for x ~ 0.23.
#   48f (x,1/8,1/8): on the -4 axis through a sodalite 4-ring (SIII).
def orbit(xyz):
    c = crystal("Li", [xyz], spacegroup=227, setting=2,
                cellpar=[A_START] * 3 + [90] * 3)
    return c.get_scaled_positions() % 1


def nearest_T(p, k=4):
    return np.sort(np.linalg.norm(mic(frac[Tidx] - p), axis=1))[:k]


# geometry sanity checks of the setting
assert np.allclose(nearest_T(np.zeros(3), 12), nearest_T(np.zeros(3), 1)[0], atol=0.05)   # D6R
assert np.allclose(nearest_T(np.full(3, 1/8), 24), nearest_T(np.full(3, 1/8), 1)[0], atol=0.05)  # sodalite

si1p = orbit((X_SI1P, X_SI1P, X_SI1P))
sii = orbit((X_SII, X_SII, X_SII))
assert len(si1p) == 32 and len(sii) == 32

# 4-ring centre on the axis from the sodalite centre: where the mean distance
# to the 4 nearest T is smallest; SIII starts D_SIII A beyond it (supercage side)
scan = np.linspace(1/8 + 3.0 / A_START, 1/8 + 5.5 / A_START, 501)
xr = scan[np.argmin([nearest_T(np.array([x, 1/8, 1/8])).mean() for x in scan])]
X_RING = float(xr)
X_SIII = X_RING + D_SIII / A_START
siii = orbit((X_SIII, 1/8, 1/8))
assert len(siii) == 48
print(f"4-ring centre x = {X_RING:.4f}, SIII start x = {X_SIII:.4f}")
# which sodalite cage each SIII site belongs to (6 per cage)
soda = orbit((1/8, 1/8, 1/8))
assert len(soda) == 8
cage_of = np.array([np.argmin(np.linalg.norm(mic(soda - p), axis=1)) for p in siii])
assert (np.bincount(cage_of) == 6).all()

# ---------------------------------------------------------------- SIII MC
fixed = np.vstack([si1p, sii])


# precompute SIII-SIII and SIII-fixed interaction matrices (periodic images
# in a 3x3x3 block with 1/r, cut at 30 A: adequate for ranking arrangements)
shifts = np.array([[i, j, k] for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1)])


def inter(Pa, Pb):
    M = np.zeros((len(Pa), len(Pb)))
    for s in shifts:
        d = (Pb[None, :, :] + s - Pa[:, None, :]) @ cell
        r = np.linalg.norm(d, axis=2)
        with np.errstate(divide="ignore"):
            w = np.where((r > 1e-6) & (r < 30.0), 1.0 / r, 0.0)
        M += w
    return M


J = inter(siii, siii)
h = inter(siii, fixed).sum(axis=1)


def energy(occ):
    o = occ.astype(float)
    return 0.5 * o @ J @ o + h @ o


found = {}
for run in range(40):
    occ = np.zeros(48, bool)
    occ[rng.choice(48, 32, replace=False)] = True
    E = energy(occ)
    T = 0.05
    for step in range(20000):
        i = rng.choice(np.where(occ)[0])
        j = rng.choice(np.where(~occ)[0])
        occ[i], occ[j] = False, True
        E2 = energy(occ)
        if E2 <= E or rng.random() < np.exp(-(E2 - E) / T):
            E = E2
        else:
            occ[i], occ[j] = True, False
        T = max(1e-4, T * 0.9995)
    key = tuple(np.where(occ)[0])
    found[key] = E

ranked = sorted(found.items(), key=lambda kv: kv[1])
print(f"{len(ranked)} distinct SIII arrangements from 40 annealing runs")
for k, E in ranked[:8]:
    print(f"  E = {E:.4f} (1/A units)")


# ---------------------------------------------------------------- write
def per_cage_counts(sel):
    # how many SIII Li per sodalite cage (6 sites per cage)
    return sorted(np.bincount(cage_of[list(sel)], minlength=8).tolist())


meta = {"x_4ring": X_RING, "x_SIII_start": X_SIII, "a_start": A_START, "x_SI'": X_SI1P, "x_SII": X_SII, "d_SIII_start_A": D_SIII,
        "candidates": []}
for n, (sel, E) in enumerate(ranked[:NCAND]):
    li = np.vstack([si1p, sii, siii[list(sel)]])
    labels = ["SI'"] * 32 + ["SII"] * 32 + ["SIII"] * 32
    at = Atoms(symbols=list(symbols) + ["Li"] * 96, cell=cell, pbc=True,
               scaled_positions=np.vstack([frac, li]))
    at.info["li_sites"] = ",".join(labels)
    write(OUT / f"LiLSX_init_c{n}.extxyz", at)
    meta["candidates"].append({"id": n, "coulomb_E": E, "siii_sites": list(map(int, sel)),
                               "siii_per_sodalite": per_cage_counts(sel)})
    print(f"c{n}: E={E:.4f}  SIII per sodalite cage {per_cage_counts(sel)}")

# min Li-framework and Li-Li distances for the first candidate
at = read(OUT / "LiLSX_init_c0.extxyz")
D = at.get_all_distances(mic=True)
li = np.where(at.numbers == 3)[0]
fw = np.where(at.numbers != 3)[0]
print("min Li-O  :", D[np.ix_(li, fw)].min().round(3))
print("min Li-Li :", (D[np.ix_(li, li)] + np.eye(96) * 99).min().round(3))
print("formula   :", at.get_chemical_formula())
(OUT / "build_meta.json").write_text(json.dumps(meta, indent=1))
np.save(OUT / "siii_sites_frac.npy", siii)
