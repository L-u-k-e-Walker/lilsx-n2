"""03_analyze.py - pick the lowest-energy relaxed Li-LSX candidate, check it
against experiment, and write the final CIFs.

Outputs (structure/):
  LiLSX_P1_MACE.cif        relaxed, fully ordered cell (672 atoms, P1 listing);
                           use this for simulations (GCMC, MD, DFT)
  LiLSX_Fd-3_averaged.cif  the same model averaged onto Fd-3 (origin 2):
                           Si, Al, O1-O4, Li SI' 32e, Li SII 32e, Li SIII 48f
                           (occupancy 2/3); crystallographic description
  results/candidates.csv   energies and cells of the relaxed candidates
  results/sites_Fd-3.csv   averaged asymmetric unit + spread over the orbit
  results/geometry.json    T-O, T-O-T, Li-O, Li coordination, comparisons
"""
import csv
import json
from pathlib import Path

import numpy as np
import spglib
from ase.io import read, write
from ase.spacegroup import Spacegroup, crystal

HERE = Path(__file__).parent
B = HERE / "build"
R = HERE / "results"
R.mkdir(exist_ok=True)
meta = json.loads((B / "build_meta.json").read_text())

# ------------------------------------------------------------ candidates
rows = []
for c in meta["candidates"]:
    f = B / f"relax_c{c['id']}.json"
    if not f.exists():
        continue
    s = json.loads(f.read_text())
    rows.append({"candidate": c["id"], "coulomb_E": round(c["coulomb_E"], 4),
                 "E_relaxed_eV": round(s["E_relaxed_eV"], 4),
                 "a": round(s["cell_lengths_A"][0], 4), "b": round(s["cell_lengths_A"][1], 4),
                 "c": round(s["cell_lengths_A"][2], 4),
                 "alpha": round(s["cell_angles_deg"][0], 3), "beta": round(s["cell_angles_deg"][1], 3),
                 "gamma": round(s["cell_angles_deg"][2], 3),
                 "V_A3": round(s["volume_A3"], 2), "max_force": round(s["max_force_eV_A"], 4),
                 "converged": s["converged"]})
emin = min(r["E_relaxed_eV"] for r in rows)
for r in rows:
    r["dE_meV_per_cell"] = round(1000 * (r["E_relaxed_eV"] - emin), 1)
with open(R / "candidates.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]))
    w.writeheader()
    w.writerows(rows)
best = min(rows, key=lambda r: r["E_relaxed_eV"])["candidate"]
print("candidates:")
for r in rows:
    print("  ", r)
print("best candidate:", best)

init = read(B / f"LiLSX_init_c{best}.extxyz")
fin = read(B / f"LiLSX_c{best}_relaxed.extxyz")
labels_li = init.info["li_sites"].split(",")
nfw = 576
assert len(fin) == 672

# ------------------------------------------------------------ geometry
D = fin.get_all_distances(mic=True)
Z = fin.numbers
iSi, iAl, iO, iLi = (np.where(Z == z)[0] for z in (14, 13, 8, 3))


def bonded(i, pool, rmax):
    return pool[D[i, pool] < rmax]


geo = {}
for name, idx in (("Si", iSi), ("Al", iAl)):
    d = np.concatenate([D[i, bonded(i, iO, 2.0)] for i in idx])
    assert len(d) == 4 * len(idx), name
    geo[f"{name}-O_mean_A"] = round(float(d.mean()), 4)
    geo[f"{name}-O_range_A"] = [round(float(d.min()), 4), round(float(d.max()), 4)]
tot = []
for o in iO:
    t = bonded(o, np.concatenate([iSi, iAl]), 2.0)
    assert len(t) == 2
    tot.append(fin.get_angle(t[0], o, t[1], mic=True))
geo["T-O-T_mean_deg"] = round(float(np.mean(tot)), 2)
geo["T-O-T_range_deg"] = [round(float(np.min(tot)), 2), round(float(np.max(tot)), 2)]

# Li: coordination (O within 2.5 A), shortest Li-O, displacement from start
li_rows = []
for k, i in enumerate(iLi):
    dO = np.sort(D[i, iO])
    v = (fin.get_scaled_positions()[i] - init.get_scaled_positions()[i])
    v -= np.round(v)
    li_rows.append({"site": labels_li[k], "cn_2.5": int((dO < 2.5).sum()),
                    "LiO_min": float(dO[0]), "LiO_mean_cn": float(dO[dO < 2.5].mean()),
                    "disp_A": float(np.linalg.norm(v @ fin.cell[:]))})
geo["Li"] = {}
for s in ("SI'", "SII", "SIII"):
    L = [r for r in li_rows if r["site"] == s]
    cn = np.bincount([r["cn_2.5"] for r in L])
    geo["Li"][s] = {"n": len(L),
                    "CN(O<2.5A)_histogram": {int(i): int(n) for i, n in enumerate(cn) if n},
                    "Li-O_shortest_mean_A": round(float(np.mean([r["LiO_min"] for r in L])), 3),
                    "Li-O_mean_within_2.5A": round(float(np.mean([r["LiO_mean_cn"] for r in L])), 3),
                    "displacement_from_start_A_mean_max": [round(float(np.mean([r["disp_A"] for r in L])), 3),
                                                           round(float(np.max([r["disp_A"] for r in L])), 3)]}
dli = D[np.ix_(iLi, iLi)] + np.eye(96) * 99
geo["Li-Li_min_A"] = round(float(dli.min()), 3)

# symmetry of the relaxed model
fw = fin[list(range(nfw))]
for tol in (1e-3, 1e-2, 5e-2, 1e-1):
    geo[f"spacegroup_framework_tol{tol}"] = spglib.get_spacegroup(
        (fw.cell[:], fw.get_scaled_positions(), fw.numbers), symprec=tol)
    geo[f"spacegroup_full_tol{tol}"] = spglib.get_spacegroup(
        (fin.cell[:], fin.get_scaled_positions(), fin.numbers), symprec=tol)

# ------------------------------------------------------------ Fd-3 average
# Orbits are defined on the (exactly Fd-3) starting structure. For every atom
# j in an orbit, find the op g with g(x_init_j) = x_init_ref, apply the same g
# to x_final_j and average; the spread is the rms deviation from the average.
sg = Spacegroup(203, setting=2)
rot, trans = sg.get_op()          # all 96 ops incl. F centring and inversion
ops = [(r, t % 1) for r, t in zip(rot, trans)]
assert len(ops) == 96
a_cub = fin.get_volume() ** (1 / 3)
xi = init.get_scaled_positions() % 1
xf = fin.get_scaled_positions()            # relaxed cell is ~cubic; frac coords used directly
cellf = fin.cell[:]


def fdist(u, v, cell):
    d = u - v
    d -= np.round(d)
    return np.linalg.norm(d @ cell)


STD_REP = {"Li SI'": np.full(3, meta["x_SI'"]), "Li SII": np.full(3, meta["x_SII"]),
           "Li SIII": np.array([meta["x_SIII_start"], 1/8, 1/8])}
site_defs = [("Si", iSi), ("Al", iAl), ("O", iO)]
li_by = {s: iLi[[k for k, l in enumerate(labels_li) if l == s]] for s in ("SI'", "SII", "SIII")}
assigned = np.zeros(len(fin), bool)
orbits = []
for name, pool in site_defs + [("Li " + s, li_by[s]) for s in li_by]:
    pool = list(pool)
    while True:
        rest = [j for j in pool if not assigned[j]]
        if not rest:
            break
        # reference member: the one in the standard representative form, so
        # the averaged coordinates read as (x,x,x) / (x,1/8,1/8) / IZA-like
        # (for SIII the standard point may be an empty 48f site, so members
        # are mapped onto the point itself rather than onto an atom)
        target = STD_REP.get(name)
        refpos = xi[rest[0]] if target is None else target % 1
        members, accum = [], []
        for j in rest:
            for r, t in ops:
                if fdist((r @ xi[j] + t) % 1, refpos, init.cell[:]) < 1e-3:
                    y = r @ xf[j] + t
                    y = refpos + ((y - refpos) - np.round(y - refpos))
                    members.append(j)
                    accum.append(y)
                    break
        for j in members:
            assigned[j] = True
        accum = np.array(accum)
        mean = accum.mean(axis=0)
        rms = float(np.sqrt((((accum - mean) @ cellf) ** 2).sum(axis=1).mean()))
        assert len(members) > 0, name
        orbits.append({"element": name.split()[0], "kind": name, "n": len(members),
                       "x": mean, "x_init": refpos, "rms_A": rms, "members": accum})
assert assigned.all()

# label O sites after the IZA FAU O1-O4 (Fd-3m) they come from
# O1-O4 coordinates read from sources/FAU_IZA.cif
txt = (HERE / "sources" / "FAU_IZA.cif").read_text().splitlines()
iza = {}
for line in txt:
    p = line.split()
    if len(p) >= 5 and p[0] in ("O1", "O2", "O3", "O4"):
        iza[p[0]] = tuple(float(v.split("(")[0]) for v in p[2:5])
iza_orb = {k: crystal("O", [v], spacegroup=227, setting=2, cellpar=[24.345] * 3 + [90] * 3)
           .get_scaled_positions() for k, v in iza.items()}
for o in orbits:
    if o["element"] == "O":
        best_lab = min(iza_orb, key=lambda k: min(fdist(o["x_init"], q, np.eye(3)) for q in iza_orb[k]))
        o["label"] = best_lab
    elif o["element"] == "Li":
        o["label"] = "Li" + o["kind"].split()[1].replace("'", "p")
    else:
        o["label"] = o["element"]
orbits.sort(key=lambda o: ["Si", "Al", "O1", "O2", "O3", "O4", "LiSIp", "LiSII", "LiSIII"].index(o["label"]))

# SI' / SII lie on the 3-fold axis (x,x,x); SIII on the -4 axis (x,1/8,1/8).
# Snap the averaged position onto its special position (the average already
# sits on it to within the rms; snapping removes the residual noise).
for o in orbits:
    x = o["x"]
    if o["label"] in ("LiSIp", "LiSII"):
        m = x.mean()
        o["x"] = np.array([m, m, m])
    if o["label"] in ("LiSIp", "LiSII"):
        # individual Li: distance from the 3-fold axis
        y = o["members"]
        u = np.ones(3) / np.sqrt(3)
        o["off_axis_each_A"] = np.linalg.norm((y - np.outer(y @ u, u)) * a_cub, axis=1)
    if o["label"] == "LiSIII":
        o["off_axis_A"] = float(np.hypot(x[1] - 0.125, x[2] - 0.125) * a_cub)
        y = o["members"]
        o["off_axis_each_A"] = np.hypot(y[:, 1] - 0.125, y[:, 2] - 0.125) * a_cub
        o["x_each"] = y[:, 0]
        o["x"] = np.array([x[0], 0.125, 0.125])
occ = {"LiSIII": 32 / 48}

# check: expanding the averaged asymmetric unit gives the right counts
for o in orbits:
    exp = crystal(o["element"], [o["x"]], spacegroup=203, setting=2,
                  cellpar=[a_cub] * 3 + [90] * 3)
    o["mult"] = len(exp)
    o["wyckoff"] = {32: "32e", 48: "48f", 96: "96g"}[len(exp)]
total = sum(o["mult"] * occ.get(o["label"], 1) for o in orbits)
counts = {}
for o in orbits:
    counts[o["element"]] = counts.get(o["element"], 0) + o["mult"] * occ.get(o["label"], 1)
assert counts == {"Si": 96, "Al": 96, "O": 384, "Li": 96}, counts

with open(R / "sites_Fd-3.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["label", "element", "wyckoff", "x", "y", "z", "occupancy", "n_atoms_averaged",
                "rms_spread_A"])
    for o in orbits:
        w.writerow([o["label"], o["element"], o["wyckoff"], *np.round(o["x"], 5),
                    round(occ.get(o["label"], 1), 4), o["n"], round(o["rms_A"], 3)])

# ------------------------------------------------------------ experiment
def cod_sites(path):
    out = {}
    for line in Path(path).read_text().splitlines():
        p = line.split()
        if len(p) == 7 and p[1].rstrip("+-0123456789") in ("Li", "O", "Si", "Al"):
            try:
                out[p[0]] = tuple(float(v) for v in p[2:5])
            except ValueError:
                pass
    return out


exp = {"Herden1982_LiX_SiAl1.23 (COD 1530083)": ("sources/1530083.cif", 24.71399),
       "Forano1989_anhydrous_LiY_SiAl2.37 (COD 1529943)": ("sources/1529943.cif", 24.4498)}
comp = {"model_a_cubic_equiv_A": round(a_cub, 4)}
# on-axis distances: SI' from D6R centre (0,0,0), SII from sodalite centre (1/8,1/8,1/8)
mod = {o["label"]: o["x"] for o in orbits}
comp["model"] = {"SI'_from_D6R_centre_A": round(np.sqrt(3) * mod["LiSIp"][0] * a_cub, 3),
                 "SII_from_sodalite_centre_A": round(np.sqrt(3) * (mod["LiSII"][0] - 0.125) * a_cub, 3),
                 "x_SI'": round(float(mod["LiSIp"][0]), 4), "x_SII": round(float(mod["LiSII"][0]), 4)}
for k, (p, a) in exp.items():
    s = cod_sites(HERE / p)
    li = sorted([v for n, v in s.items() if n.startswith("Li")], key=lambda v: v[0])
    comp[k] = {"a_A": a, "x_SI'": li[0][0], "x_SII": li[1][0],
               "SI'_from_D6R_centre_A": round(np.sqrt(3) * li[0][0] * a, 3),
               "SII_from_sodalite_centre_A": round(np.sqrt(3) * (li[1][0] - 0.125) * a, 3)}
    # framework O: distance (A, after scaling to the model cell) between each
    # model O site and the nearest experimental O in the Fd-3m orbit
    dev = {}
    for o in orbits:
        if o["element"] != "O":
            continue
        # COD labels O1-O4 do not follow the IZA numbering: match against all
        best_d, best_l = 1e9, None
        for lab, v in s.items():
            if not lab.startswith("O"):
                continue
            e = crystal("O", [v], spacegroup=227, setting=2,
                        cellpar=[a] * 3 + [90] * 3).get_scaled_positions()
            d = min(fdist(o["x"], q, np.eye(3) * a_cub) for q in e)
            if d < best_d:
                best_d, best_l = d, lab
        dev[o["label"]] = {"matched_exp_label": best_l, "offset_A": round(best_d, 3)}
    comp[k]["O_site_offset_from_model_A"] = dev
geo["comparison"] = comp
o3 = [o for o in orbits if o["label"] == "LiSIII"][0]
geo["SIII_average_off_axis_before_snapping_A"] = round(o3["off_axis_A"], 3)
geo["Li_individual_off_axis_A_mean_max"] = {
    o["label"]: [round(float(o["off_axis_each_A"].mean()), 3), round(float(o["off_axis_each_A"].max()), 3)]
    for o in orbits if "off_axis_each_A" in o}
geo["SIII_individual_x_range"] = [round(float(o3["x_each"].min()), 4), round(float(o3["x_each"].max()), 4)]
geo["SIII_4ring_centre_x"] = round(meta["x_4ring"], 4)
geo["SIII_beyond_4ring_centre_A"] = round(float((o3["x"][0] - meta["x_4ring"]) * a_cub), 3)
rc = json.loads((B / f"relax_c{best}.json").read_text())
geo["cubic_cell_pressure_GPa"] = round(-float(np.mean(rc["cubic"]["stress_GPa"][:3])), 4)
geo["cubic_cell_max_abs_stress_component_GPa"] = round(float(np.max(np.abs(rc["cubic"]["stress_GPa"]))), 4)
geo["free_cell"] = {k: rc["free_cell"][k] for k in ("cell_lengths_A", "cell_angles_deg", "volume_A3")}
geo["free_minus_cubic_meV_per_cell"] = round(1000 * (rc["free_cell"]["E_eV"] - rc["cubic"]["E_eV"]), 1)
geo["orbit_rms_spread_A"] = {o["label"]: round(o["rms_A"], 3) for o in orbits}
geo["best_candidate"] = best
(R / "geometry.json").write_text(json.dumps(geo, indent=1))
print(json.dumps(geo, indent=1))

# ------------------------------------------------------------ CIFs
HEAD = """# Li-LSX model: Li96 [Al96 Si96 O384], FAU, Si/Al = 1 (strict alternation)
# Built from the IZA-SC FAU framework (Baerlocher & McCusker), Li at SI' (32),
# SII (32), SIII (32 of 48) following dehydrated Li-LSX literature, relaxed with
# MACE-MP-0 medium + D3(BJ) (float32), cubic cell. NOT an experimental structure.
# See lilsx-n2/structure/README.md for sources, method and validation.
"""
p1 = HERE / "LiLSX_P1_MACE.cif"
# labels = element only (RASPA reads labels as pseudo-atom types); the Li
# crystallographic site of every atom is in results/li_sites_P1.csv
write(p1, fin, format="cif", labels=[fin.get_chemical_symbols()])
p1.write_text(HEAD + "# Ordered model, all 672 atoms listed (P1). Cell relaxed under a cubic\n"
              "# constraint (hydrostatic strain only, pressure ~0); atom order: framework\n"
              "# 0-575, then Li SI' (32), SII (32), SIII (32).\n" + p1.read_text())
with open(R / "li_sites_P1.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["atom_index_0based", "site", "x", "y", "z", "n_O_within_2.5A", "Li-O_min_A"])
    fp = fin.get_scaled_positions()
    for k, i in enumerate(iLi):
        w.writerow([int(i), labels_li[k], *np.round(fp[i], 5), li_rows[k]["cn_2.5"],
                    round(li_rows[k]["LiO_min"], 3)])

lines = [HEAD + "# Fd-3 average of the ordered model: SIII is 48f with occupancy 2/3.\n",
         "data_LiLSX_Fd-3_averaged",
         f"_chemical_formula_sum 'Li96 Al96 Si96 O384'",
         f"_cell_length_a {a_cub:.4f}", f"_cell_length_b {a_cub:.4f}", f"_cell_length_c {a_cub:.4f}",
         "_cell_angle_alpha 90", "_cell_angle_beta 90", "_cell_angle_gamma 90",
         "_symmetry_space_group_name_H-M 'F d -3'", "_symmetry_Int_Tables_number 203",
         "_space_group.IT_coordinate_system_code '2'", "", "loop_", "_symmetry_equiv_pos_as_xyz"]


def op_str(r, t):
    out = []
    for i in range(3):
        s = ""
        for j, v in enumerate("xyz"):
            if r[i, j] == 1:
                s += "+" + v
            elif r[i, j] == -1:
                s += "-" + v
        tt = t[i] % 1
        if abs(tt) > 1e-6:
            from fractions import Fraction
            s += "+" + str(Fraction(tt).limit_denominator(8))
        out.append(s.lstrip("+"))
    return "'" + ",".join(out) + "'"


seen = set()
for r, t in ops:
    s = op_str(r, t)
    if s not in seen:
        seen.add(s)
        lines.append(s)
assert len(seen) == 96, len(seen)
lines += ["", "loop_", "_atom_site_label", "_atom_site_type_symbol", "_atom_site_Wyckoff_symbol",
          "_atom_site_fract_x", "_atom_site_fract_y", "_atom_site_fract_z", "_atom_site_occupancy"]
for o in orbits:
    lines.append(f"{o['label']:7s} {o['element']:2s} {o['wyckoff']} " +
                 " ".join(f"{v:.5f}" for v in o["x"]) + f" {occ.get(o['label'], 1):.6f}")
(HERE / "LiLSX_Fd-3_averaged.cif").write_text("\n".join(lines) + "\n")

print("wrote", p1.name, "and LiLSX_Fd-3_averaged.cif")

# ------------------------------------------------------------ 6-ring heights
# signed distance of SI' / SII Li from the plane of its 3 nearest O
# (positive = into the sodalite cage for SI', into the supercage for SII)
pos = fin.get_positions()
heights = {}
for s, ref_centre in (("SI'", np.full(3, 1/8)), ("SII", np.full(3, 3/8))):
    h = []
    for k, i in enumerate(iLi):
        if labels_li[k] != s:
            continue
        o3 = iO[np.argsort(D[i, iO])[:3]]
        v = [fin.get_distance(i, o, mic=True, vector=True) for o in o3]
        P = pos[i] + np.array(v)                  # the 3 O, unwrapped around Li
        n = np.cross(P[1] - P[0], P[2] - P[0])
        n /= np.linalg.norm(n)
        d = float((pos[i] - P.mean(axis=0)) @ n)
        # orient: towards the nearest sodalite centre (SI') / away from it (SII)
        soda = crystal("X", [(1/8, 1/8, 1/8)], spacegroup=227, setting=2,
                       cellpar=[1, 1, 1, 90, 90, 90]).get_scaled_positions()
        fi = fin.get_scaled_positions()[i]
        dd = soda - fi
        dd -= np.round(dd)
        c = (dd[np.argmin(np.linalg.norm(dd @ fin.cell[:], axis=1))] @ fin.cell[:])
        to_soda = np.sign(n @ c)
        h.append(d * to_soda * (1 if s == "SI'" else -1))
    heights[s] = {"mean_A": round(float(np.mean(h)), 3), "min_A": round(float(np.min(h)), 3),
                  "max_A": round(float(np.max(h)), 3)}
geo["Li_height_from_6ring_O_plane_A"] = heights


def exp_heights(path, a):
    """same quantity for an experimental Fd-3m structure (COD)"""
    st = cod_sites(HERE / path)
    O = np.vstack([crystal("O", [v], spacegroup=227, setting=2, cellpar=[a] * 3 + [90] * 3)
                   .get_scaled_positions() for n, v in st.items() if n.startswith("O")])
    li = sorted([v for n, v in st.items() if n.startswith("Li")], key=lambda v: v[0])
    out = {}
    for s_, x, sodac, sgn in (("SI'", li[0], np.full(3, 1/8), 1), ("SII", li[1], np.full(3, 1/8), -1)):
        x = np.array(x)
        d = O - x
        d -= np.round(d)
        dc = d @ (np.eye(3) * a)
        P = dc[np.argsort(np.linalg.norm(dc, axis=1))[:3]]
        n = np.cross(P[1] - P[0], P[2] - P[0])
        n /= np.linalg.norm(n)
        h = float((0 - P.mean(axis=0)) @ n)
        to_soda = np.sign(n @ ((sodac - x) * a))
        out[s_] = {"height_A": round(h * to_soda * sgn, 3),
                   "Li-O_3nearest_A": [round(float(v), 3) for v in np.sort(np.linalg.norm(dc, axis=1))[:3]]}
    return out


for k, (p, a) in exp.items():
    geo["comparison"][k]["Li_height_from_6ring_O_plane_A"] = exp_heights(p, a)
    print(k, geo["comparison"][k]["Li_height_from_6ring_O_plane_A"])
print("6-ring heights:", heights)

# ------------------------------------------------------------ CIF round trip
from pymatgen.core import Structure
p1s = Structure.from_file(str(p1))
avg = Structure.from_file(str(HERE / "LiLSX_Fd-3_averaged.cif"))
geo["cif_check"] = {"P1_atoms": len(p1s), "P1_formula": p1s.composition.formula,
                    "Fd-3_sites_expanded": len(avg),
                    "Fd-3_composition": avg.composition.formula,
                    "Fd-3_spacegroup_pymatgen": avg.get_space_group_info()[0] if avg.is_ordered else "disordered (SIII 2/3)"}
print("cif check:", geo["cif_check"])
(R / "geometry.json").write_text(json.dumps(geo, indent=1))
