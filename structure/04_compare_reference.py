"""04_compare_reference.py - compare the Li-LSX model with a published structure.

Reference: ../reference_structure/1540962.cif (COD 1540962)
  Shepelev, Anderson, Smolin, Zeolites 10 (1990) 61-63, "Crystal structure of a
  partially lithium-exchanged X zeolite in hydrated (25 C) and dehydrated
  (275 C) states". Li41.28 Na24.32 [Al92 Si100 O384] (H2O)8.64, Fd-3 origin 2,
  a = 24.75 A, Si/Al ordered (Si1 full Si; Al1 = 0.958 Al / 0.042 Si).

Model: results/sites_Fd-3.csv (Fd-3 average of LiLSX_P1_MACE.cif, a = 24.7638 A).

Both are Fd-3 origin 2, but Fd-3 sits in Fd-3m with index 2, so the same
structure can be written with Si and Al swapped between the two T orbits (plus
other equivalent origin/axis choices). The model is therefore first mapped
with the Fd-3m operation that best overlays its framework on the reference.

Outputs: results/compare_shepelev1990.json, results/compare_shepelev1990_sites.csv
"""
import csv
import json
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.spacegroup import Spacegroup

HERE = Path(__file__).parent
REF = HERE.parent / "reference_structure" / "1540962.cif"
R = HERE / "results"


# ------------------------------------------------------------ parsing
def parse_ops(lines):
    ops = []
    for s in lines:
        comps = s.strip().strip("'").split(",")
        rot = np.zeros((3, 3))
        tr = np.zeros(3)
        for i, c in enumerate(comps):
            env0 = {"x": 0.0, "y": 0.0, "z": 0.0}
            tr[i] = eval(c, {}, env0)
            for j, v in enumerate("xyz"):
                env = dict(env0)
                env[v] = 1.0
                rot[i, j] = eval(c, {}, env) - tr[i]
        ops.append((np.rint(rot).astype(int), tr % 1))
    return ops


def read_cif(path):
    txt = Path(path).read_text().splitlines()
    i = txt.index("_symmetry_equiv_pos_as_xyz") + 1
    op_lines = []
    while not txt[i].startswith(("loop_", "_")):
        op_lines.append(txt[i])
        i += 1
    sites = {}
    for line in txt:
        p = line.split()
        if len(p) == 7 and p[1].rstrip("+-0123456789") in ("Li", "Na", "O", "Si", "Al"):
            sites[p[0]] = {"el": p[1].rstrip("+-0123456789"), "x": np.array(list(map(float, p[2:5]))),
                           "occ": float(p[5])}
    a = float([l for l in txt if l.startswith("_cell_length_a")][0].split()[1])
    return a, parse_ops(op_lines), sites


a_ref, ops_ref, ref = read_cif(REF)
assert len(ops_ref) == 96
sg203 = Spacegroup(203, setting=2)
ops_model = [(r, t % 1) for r, t in zip(*sg203.get_op())]


def opset(ops):
    return {(tuple(r.ravel()), tuple(np.round(t % 1, 4) % 1)) for r, t in ops}


assert opset(ops_ref) == opset(ops_model), "reference is not in the same Fd-3 origin-2 setting"

model = {}
with open(R / "sites_Fd-3.csv") as fh:
    for row in csv.DictReader(fh):
        model[row["label"]] = {"el": row["element"], "x": np.array([float(row[k]) for k in "xyz"]),
                               "occ": float(row["occupancy"]), "wyck": row["wyckoff"]}
geo = json.loads((R / "geometry.json").read_text())
a_mod = geo["comparison"]["model_a_cubic_equiv_A"]
a_avg = 0.5 * (a_ref + a_mod)


def orbit(x, ops=ops_ref):
    pts = np.array([(r @ x + t) % 1 for r, t in ops])
    keep = []
    for p in pts:
        if not any(np.linalg.norm(((p - q + 0.5) % 1) - 0.5) < 1e-4 for q in keep):
            keep.append(p)
    return np.array(keep)


def fd(u, v):
    d = u - v
    return d - np.round(d)


def min_dist(x, orb, a=a_avg):
    d = fd(orb, x)
    k = np.argmin(np.linalg.norm(d, axis=1))
    return float(np.linalg.norm(d[k]) * a), orb[k]


ref_orb = {k: orbit(v["x"]) for k, v in ref.items()}

# ------------------------------------------------------------ best overlay
# T sites of the reference: Si1 = Si, Al1 (0.958 Al) = Al
ref_T = {"Si": "Si1", "Al": "Al1"}
ref_O = [k for k in ref if ref[k]["el"] == "O" and k != "O5"]   # O5 = water at 16c
fw_model = ["Si", "Al", "O1", "O2", "O3", "O4"]
sg227 = Spacegroup(227, setting=2)
best = None
for r, t in zip(*sg227.get_op()):
    score = 0.0
    for lab in fw_model:
        y = (r @ model[lab]["x"] + t) % 1
        if lab in ref_T:
            d, _ = min_dist(y, ref_orb[ref_T[lab]])
        else:
            d = min(min_dist(y, ref_orb[o])[0] for o in ref_O)
        score += d
    if best is None or score < best[0]:
        best = (score, r, t)
_, Rg, Tg = best
mapped = {k: {**v, "x": (Rg @ v["x"] + Tg) % 1} for k, v in model.items()}

# ------------------------------------------------------------ site table
rows = []
used = set()
for lab, m in mapped.items():
    if lab in ref_T:
        cands = [ref_T[lab]]
    elif m["el"] == "O":
        cands = ref_O
    else:   # Li sites: compare with every reference cation site
        cands = [k for k in ref if ref[k]["el"] in ("Li", "Na")]
    dists = sorted((min_dist(m["x"], ref_orb[c])[0], c) for c in cands)
    d, c = dists[0]
    _, img = min_dist(m["x"], ref_orb[c])
    rows.append({"model_site": lab, "model_xyz": np.round(m["x"], 4).tolist(),
                 "ref_site": c, "ref_el": ref[c]["el"], "ref_occ": ref[c]["occ"],
                 "ref_xyz_image": np.round(img, 4).tolist(), "offset_A": round(d, 3)})

# ------------------------------------------------------------ geometry
def expand(sites, a, labels):
    sym, pos, lab = [], [], []
    for k in labels:
        for p in orbit(sites[k]["x"]):
            sym.append(sites[k]["el"])
            pos.append(p)
            lab.append(k)
    at = Atoms(sym, scaled_positions=pos, cell=np.eye(3) * a, pbc=True)
    return at, np.array(lab)


def framework_geometry(sites, a, Tlab, Olab):
    at, lab = expand(sites, a, list(Tlab) + list(Olab))
    D = at.get_all_distances(mic=True)
    out = {}
    for k in Tlab:
        i = np.where(lab == k)[0][0]
        d = np.sort(D[i][np.isin(lab, Olab)])[:4]
        out[f"{k}-O_mean_A"] = round(float(d.mean()), 4)
        out[f"{k}-O_each_A"] = np.round(d, 3).tolist()
    isT = np.isin(lab, Tlab)
    for k in Olab:
        i = np.where(lab == k)[0][0]
        t = np.where(isT)[0][np.argsort(D[i][isT])[:2]]
        out[f"T-{k}-T_deg"] = round(float(at.get_angle(t[0], i, t[1], mic=True)), 2)
    return out


def cation_geometry(sites, a, cat, Olab):
    """nearest framework O distances and height above the 6-ring O plane
    (SI'/SII-type sites on the 3-fold axis) for one cation site"""
    at, lab = expand(sites, a, list(Olab))
    x = sites[cat]["x"]
    d = fd(at.get_scaled_positions(), x) @ at.cell[:]
    r = np.linalg.norm(d, axis=1)
    k = np.argsort(r)[:4]
    res = {"nearest_O_A": np.round(r[k], 3).tolist(), "nearest_O_sites": lab[k].tolist()}
    if np.allclose(x, x[0], atol=1e-6):            # on the 3-fold axis
        P = d[k[:3]]
        n = np.cross(P[1] - P[0], P[2] - P[0])
        n /= np.linalg.norm(n)
        h = float(-(P.mean(axis=0)) @ n)
        # orient towards the sodalite centre (1/8,1/8,1/8) for x < 1/8 (SI'),
        # away from it for x > 1/8 (SII): positive = away from the 6-ring into
        # the cage the site belongs to
        c = (np.full(3, 1/8) - x) * a
        s = np.sign(n @ c) * (1 if x[0] < 1/8 else -1)
        res["height_from_6ring_O_plane_A"] = round(h * s, 3)
    return res


Olab_ref = ref_O
geo_ref = framework_geometry(ref, a_ref, ["Si1", "Al1"], Olab_ref)
geo_mod = framework_geometry(mapped, a_mod, ["Si", "Al"], ["O1", "O2", "O3", "O4"])
cat_ref = {k: cation_geometry(ref, a_ref, k, Olab_ref) for k in ref if ref[k]["el"] in ("Li", "Na")}
# geometry is invariant under the overlay op; use the unmapped (x,x,x)-form sites
cat_mod = {k: cation_geometry(model, a_mod, k, ["O1", "O2", "O3", "O4"]) for k in ("LiSIp", "LiSII", "LiSIII")}

# T-O-T per O site: pair model O label with the reference O it overlays
o_pair = {r_["model_site"]: r_["ref_site"] for r_ in rows if r_["model_site"].startswith("O")}
tot = [{"model_O": mo, "ref_O": ro, "model_TOT": geo_mod[f"T-{mo}-T_deg"], "ref_TOT": geo_ref[f"T-{ro}-T_deg"]}
       for mo, ro in o_pair.items()]

out = {
    "reference": {"cod": 1540962, "citation": "Shepelev, Anderson, Smolin, Zeolites 10 (1990) 61-63",
                  "composition": "Li41.28 Na24.32 [Al92 Si100 O384] (H2O)8.64", "Si/Al": round(100 / 92, 3),
                  "a_A": a_ref, "space_group": "Fd-3 (origin 2)",
                  "cations_located_per_cell": round(sum(ref[k]["occ"] * len(ref_orb[k]) for k in ref
                                                       if ref[k]["el"] in ("Li", "Na")), 2),
                  "cations_needed_per_cell": 92,
                  "site_multiplicity_and_count": {k: [len(ref_orb[k]), round(ref[k]["occ"] * len(ref_orb[k]), 2)]
                                                  for k in ref}},
    "model": {"a_A": a_mod, "composition": "Li96 [Al96 Si96 O384]"},
    "a_difference_A": round(a_mod - a_ref, 4),
    "a_difference_pct": round(100 * (a_mod - a_ref) / a_ref, 3),
    "overlay_op": {"rotation": Rg.tolist(), "translation": np.round(Tg, 4).tolist(),
                   "framework_offset_sum_A": round(best[0], 3)},
    "sites": rows,
    "framework_geometry": {"reference": geo_ref, "model": geo_mod},
    "T-O-T_by_site": tot,
    "cation_geometry": {"reference": cat_ref, "model": cat_mod},
}
(R / "compare_shepelev1990.json").write_text(json.dumps(out, indent=1))
with open(R / "compare_shepelev1990_sites.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["model_site", "model_x", "model_y", "model_z", "ref_site", "ref_element",
                "ref_occupancy", "ref_x", "ref_y", "ref_z", "offset_A"])
    for r_ in rows:
        w.writerow([r_["model_site"], *r_["model_xyz"], r_["ref_site"], r_["ref_el"], r_["ref_occ"],
                    *r_["ref_xyz_image"], r_["offset_A"]])
print(json.dumps(out, indent=1))
