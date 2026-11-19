"""make_trajectory.py - scripted N2 adsorption storyboard (env viz). Not MD.

120 frames at 15 fps (t = 0 ... 7.93 s), framework static:
  0.0-2.5 s  N2 drifts and tumbles in from near the supercage centre
  2.5-3.5 s  rotates end-on and docks at the bound geometry (eased)
  3.5-6.0 s  bound: small smooth wobble of the Li-side N and libration
  6.0-8.0 s  lifts off end-on and drifts away, back towards the cage centre

Seamless loop: every motion term has period 8 s (tumble rates and noise
frequencies are integer multiples of 1/8 Hz, and the path ends where it
started with the same velocity), so the GIF repeats without a jump.

Smooth noise only: sums of low-frequency sinusoids with fixed random phases,
never per-frame random numbers.

Geometry:
  bound pose from ../binding/complex_relaxed.extxyz (MACE-MP-0 + D3, N2 only):
  Li 662 (SIII), Li-N 2.063 A, Li-N-N 171.2 deg; N-N kept at the bound value.
  The molecule is placed through a reference point P_s = N_in + s*d*u, with
  s = 0.5 (centre of mass) in free flight and s = 0 (pivot on the Li-side N)
  when bound, blended with the same weight w(t) as the orientation.

Outputs: anim.extxyz (674 atoms/frame; per-frame info: time, Li-N, angle,
glow flag), trajectory_qc.csv, qc_summary.json. The script fails if any frame
breaks the distance limits below.
"""
import csv
import json
from pathlib import Path

import numpy as np
from ase.io import read, write

HERE = Path(__file__).parent
ROOT = HERE.parents[1]
NFRAMES, FPS = 120, 15
T = NFRAMES / FPS                     # 8 s period
t_all = np.arange(NFRAMES) / FPS
rng = np.random.default_rng(20261118)

# distance limits checked every frame (A)
MIN_N_FRAMEWORK = 2.9                 # N to any framework O/Si/Al
MIN_N_OTHER_LI = 2.9                  # N to any Li except the active one
MIN_N_ACTIVE_LI = 1.95                # Li-side N to the active Li
GLOW_CUTOFF = 2.8                     # Li...N connector shown below this

# ---------------------------------------------------------------- bound pose
cx = read(ROOT / "animation" / "binding" / "complex_relaxed.extxyz")
summ = json.loads((ROOT / "animation" / "binding" / "binding_summary.json").read_text())
LI = summ["Li_index"]
NFW = len(cx) - 2
fw = cx[:NFW]
A = fw.cell[0, 0]
CENTRE = np.full(3, 3 / 8) * A
VIEW = np.ones(3) / np.sqrt(3)        # Phase-4 camera direction (camera -> centre)


def mic(v):
    return v - np.round(v / A) * A


li = fw.positions[LI]
na, nb = cx.positions[NFW], cx.positions[NFW + 1]
if np.linalg.norm(mic(nb - li)) < np.linalg.norm(mic(na - li)):
    na, nb = nb, na
n_in_b = li + mic(na - li)
n_out_b = n_in_b + mic(nb - na)
D_NN = float(np.linalg.norm(n_out_b - n_in_b))
u_b = (n_out_b - n_in_b) / D_NN
v_hat = mic(n_in_b + 0.5 * D_NN * u_b - li)
v_hat /= np.linalg.norm(v_hat)        # Li -> bound N2 centre direction


# ---------------------------------------------------------------- helpers
def smootherstep(x):
    x = np.clip(x, 0.0, 1.0)
    return x * x * x * (x * (6 * x - 15) + 10)


def qmul(a, b):
    w1, x1, y1, z1 = a
    w2, x2, y2, z2 = b
    return np.array([w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2,
                     w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2,
                     w1 * y2 - x1 * z2 + y1 * w2 + z1 * x2,
                     w1 * z2 + x1 * y2 - y1 * x2 + z1 * w2])


def qaxis(axis, ang):
    axis = np.asarray(axis, float) / np.linalg.norm(axis)
    return np.concatenate([[np.cos(ang / 2)], np.sin(ang / 2) * axis])


def qrot(q, v):
    p = np.concatenate([[0.0], v])
    qc = q * np.array([1, -1, -1, -1])
    return qmul(qmul(q, p), qc)[1:]


def qfrom_to(a, b):
    a, b = a / np.linalg.norm(a), b / np.linalg.norm(b)
    c = np.cross(a, b)
    if np.linalg.norm(c) < 1e-9:
        return np.array([1.0, 0, 0, 0])
    return qaxis(c, np.arccos(np.clip(a @ b, -1, 1)))


def slerp(q0, q1, w):
    if q0 @ q1 < 0:
        q1 = -q1
    d = np.clip(q0 @ q1, -1, 1)
    th = np.arccos(d)
    if th < 1e-6:
        return q0
    return (np.sin((1 - w) * th) * q0 + np.sin(w * th) * q1) / np.sin(th)


def smooth_noise(n_terms=3, kmax=4, dim=3):
    """periodic (8 s) low-frequency noise, unit-ish amplitude per component"""
    ks = rng.integers(1, kmax + 1, size=(dim, n_terms))
    ph = rng.uniform(0, 2 * np.pi, size=(dim, n_terms))
    amp = rng.uniform(0.5, 1.0, size=(dim, n_terms))
    norm = amp.sum(axis=1)

    def f(t):
        return np.array([(amp[i] * np.sin(2 * np.pi * ks[i] * t / T + ph[i])).sum() / norm[i]
                         for i in range(dim)])
    return f


def perp_basis(u):
    a = np.cross(u, [0, 0, 1.0])
    if np.linalg.norm(a) < 1e-3:
        a = np.cross(u, [0, 1.0, 0])
    a /= np.linalg.norm(a)
    return a, np.cross(u, a)


def hermite(keys, t):
    """C1 cubic Hermite through (time, position, velocity) keyframes"""
    for (t0, p0, v0), (t1, p1, v1) in zip(keys[:-1], keys[1:]):
        if t0 <= t <= t1:
            h = t1 - t0
            x = (t - t0) / h
            h00, h10 = 2 * x**3 - 3 * x**2 + 1, x**3 - 2 * x**2 + x
            h01, h11 = -2 * x**3 + 3 * x**2, x**3 - x**2
            return h00 * p0 + h10 * h * v0 + h01 * p1 + h11 * h * v1
    raise ValueError(t)


# ---------------------------------------------------------------- storyboard
# screen basis of the Phase-4 camera (up = z projected)
up = np.array([0, 0, 1.0]) - VIEW[2] * VIEW
up /= np.linalg.norm(up)
right = np.cross(VIEW, up)
lat = mic(li - CENTRE)
lat -= (lat @ VIEW) * VIEW            # Li offset from the centre, in the screen plane

# start/end point: just across the centre from the Li, slightly towards the camera
P0 = CENTRE - 0.15 * lat - 1.0 * VIEW + 0.6 * up
V0 = 0.6 * (lat / np.linalg.norm(lat))                        # gentle drift towards the Li side
P_PRE = li + 4.4 * v_hat                                       # pre-dock point
V_PRE = -1.4 * v_hat + 0.4 * np.cross(v_hat, VIEW)             # arriving, slightly curved
# lift-off direction: away from the active Li, tilted away from the nearest
# other Li (an SII Li 3.5 A from the outer N in the bound pose), so the
# un-docking rotation cannot swing an N into it
com_b = n_in_b + 0.5 * D_NN * u_b
others = [j for j in np.where(np.array(fw.get_chemical_symbols()) == "Li")[0] if j != LI]
j_near = min(others, key=lambda j: np.linalg.norm(mic(fw.positions[j] - com_b)))
away = mic(com_b - fw.positions[j_near])
away -= (away @ v_hat) * v_hat
d_off = v_hat + 0.7 * away / np.linalg.norm(away)
d_off /= np.linalg.norm(d_off)
P_OFF = n_in_b + 3.0 * d_off                                   # after lift-off
# leave P_OFF heading for P0 at the average speed needed to arrive by T,
# blended with the lift-off direction so the path stays a smooth curve
to_p0 = P0 - P_OFF
speed = np.linalg.norm(to_p0) / (T - 6.9)
V_OFF = speed * (0.5 * d_off + 0.5 * to_p0 / np.linalg.norm(to_p0))
keys = [(0.0, P0, V0), (2.5, P_PRE, V_PRE), (3.5, n_in_b, np.zeros(3)),
        (6.0, n_in_b, np.zeros(3)), (6.9, P_OFF, V_OFF), (T, P0, V0)]


def w_bind(t):
    """0 = free tumbling, 1 = bound orientation/pivot"""
    if t < 1.6:
        return 0.0
    if t < 3.5:
        return float(smootherstep((t - 1.6) / 1.9))
    if t < 6.3:
        return 1.0
    if t < 7.9:
        return float(1 - smootherstep((t - 6.3) / 1.6))
    return 0.0


# free tumble: rotations about two fixed axes at 2 and 1 turns per 8 s
# The tumble axes/phase are arbitrary; they are chosen by a small search
# (below) so the free orientation is close to the bound one around docking
# and undocking, which keeps the blended rotation slow.
n1 = rng.normal(size=3); n1 /= np.linalg.norm(n1)
n2 = rng.normal(size=3); n2 /= np.linalg.norm(n2)
q_ref = qfrom_to(np.array([0, 0, 1.0]), -v_hat)
q_bound = qfrom_to(np.array([0, 0, 1.0]), u_b)
tumble_wobble = smooth_noise(n_terms=2, kmax=3, dim=2)


def q_tumble(t):
    wob = tumble_wobble(t)
    a1 = 2 * np.pi * 2 * t / T + 0.35 * wob[0]
    a2 = 2 * np.pi * 1 * t / T + 0.35 * wob[1]
    return qmul(qaxis(n1, a1), qmul(qaxis(n2, a2), q_ref))


# bound libration (pivot on the Li-side N) and anchor wobble
e1, e2 = perp_basis(u_b)
lib_noise = smooth_noise(n_terms=3, kmax=6, dim=2)
LIB_AMP = np.radians(7.0)
anchor_noise = smooth_noise(n_terms=3, kmax=6, dim=3)
ANCHOR_AMP = 0.06                     # A
free_noise = smooth_noise(n_terms=3, kmax=3, dim=3)
FREE_AMP = 0.30                       # A


def frame(t):
    w = w_bind(t)
    lib = lib_noise(t) * LIB_AMP
    q_b = qmul(qaxis(e1, lib[0]), qmul(qaxis(e2, lib[1]), q_bound))
    q = slerp(q_tumble(t), q_b, w)
    u = qrot(q, np.array([0, 0, 1.0]))
    s = 0.5 * (1 - w)
    ref = hermite(keys, t) + (1 - w) * FREE_AMP * free_noise(t) + w * ANCHOR_AMP * anchor_noise(t)
    n_in = ref - s * D_NN * u
    return n_in, n_in + D_NN * u, w


# ---------------------------------------------------------------- tumble search
MIN_TUMBLE_IN, MIN_TUMBLE_OUT = 150.0, 90.0   # deg of rotation in 0-2 s and 6.8-8 s


def max_rot_and_ok():
    worst, prev_u, cin, cout = 0.0, None, 0.0, 0.0
    for t in np.append(t_all, T):
        n_in, n_out, _ = frame(t)
        u = (n_out - n_in) / D_NN
        if prev_u is not None:
            step = np.degrees(np.arccos(np.clip(u @ prev_u, -1, 1)))
            worst = max(worst, step)
            cin += step if t <= 2.0 else 0.0
            cout += step if t > 6.8 else 0.0
        prev_u = u
        N = np.array([n_in, n_out])
        d = np.linalg.norm(mic(fw.positions[None, :, :] - N[:, None, :]), axis=2)
        dl = d[:, LI].min()
        d[:, LI] = 99
        is_l = np.array(fw.get_chemical_symbols()) == "Li"
        if d[:, ~is_l].min() < MIN_N_FRAMEWORK or d[:, is_l].min() < MIN_N_OTHER_LI or dl < MIN_N_ACTIVE_LI:
            return worst, False
    return worst, (cin >= MIN_TUMBLE_IN and cout >= MIN_TUMBLE_OUT)


search_rng = np.random.default_rng(7)
best = None
for trial in range(300):
    c1 = search_rng.normal(size=3); c1 /= np.linalg.norm(c1)
    c2 = search_rng.normal(size=3); c2 /= np.linalg.norm(c2)
    cq = search_rng.normal(size=4); cq /= np.linalg.norm(cq)
    n1, n2, q_ref = c1, c2, cq
    worst, ok = max_rot_and_ok()
    if ok and (best is None or worst < best[0]):
        best = (worst, c1, c2, cq, trial)
TUMBLE_SEARCH = {"trials": 300, "min_rotation_deg_0-2s": MIN_TUMBLE_IN, "min_rotation_deg_6.8-8s": MIN_TUMBLE_OUT, "best_trial": best[4], "max_rot_deg_per_frame": round(float(best[0]), 2)}
_, n1, n2, q_ref, _ = best
print("tumble search:", TUMBLE_SEARCH)

# ---------------------------------------------------------------- build + QC
fw_pos = fw.positions
is_li = np.array(fw.get_chemical_symbols()) == "Li"
other_li = is_li.copy()
other_li[LI] = False
non_li = ~is_li
images, rows = [], []
prev = None
for i, t in enumerate(t_all):
    n_in, n_out, w = frame(t)
    N = np.array([n_in, n_out])
    d = np.linalg.norm(mic(fw_pos[None, :, :] - N[:, None, :]), axis=2)
    d_li = d[:, LI]
    k_in = int(np.argmin(d_li))
    li_n = float(d_li.min())
    near_n = N[k_in]
    far_n = N[1 - k_in]
    a = mic(li - near_n)
    b = far_n - near_n
    ang = float(np.degrees(np.arccos(np.clip(a @ b / np.linalg.norm(a) / np.linalg.norm(b), -1, 1))))
    row = {"frame": i, "time_s": round(float(t), 4), "w_bind": round(w, 4),
           "Li-N_A": round(li_n, 3), "Li-N-N_deg": round(ang, 1),
           "min_N_framework_A": round(float(d[:, non_li].min()), 3),
           "min_N_other_Li_A": round(float(d[:, other_li].min()), 3),
           "glow": int(li_n < GLOW_CUTOFF),
           "com_x": round(float((n_in[0] + n_out[0]) / 2), 3),
           "com_y": round(float((n_in[1] + n_out[1]) / 2), 3),
           "com_z": round(float((n_in[2] + n_out[2]) / 2), 3)}
    if prev is not None:
        row["step_N_max_A"] = round(float(np.linalg.norm(N - prev[0], axis=1).max()), 3)
        uu = (n_out - n_in) / D_NN
        row["step_rot_deg"] = round(float(np.degrees(np.arccos(np.clip(uu @ prev[1], -1, 1)))), 2)
    prev = (N, (n_out - n_in) / D_NN)
    rows.append(row)
    img = fw.copy()
    img.extend(cx[NFW:NFW + 2])
    img.positions[NFW] = n_in
    img.positions[NFW + 1] = n_out
    img.info = {"time_s": float(t), "frame": i, "li_n_A": li_n, "li_n_n_deg": ang,
                "glow": int(li_n < GLOW_CUTOFF), "active_li": LI, "n_inner": NFW + k_in, "w_bind": w}
    img.calc = None
    images.append(img)

# loop seam: frame 0 vs the frame that would follow frame 119
n_in_T, n_out_T, _ = frame(T)
seam = float(np.linalg.norm(np.array([n_in_T, n_out_T]) - np.array(frame(0.0)[:2]), axis=1).max())
rows[0]["step_N_max_A"] = round(float(np.linalg.norm(np.array(frame(0.0)[:2]) - np.array(frame(T - 1 / FPS)[:2]), axis=1).max()), 3)

with open(HERE / "trajectory_qc.csv", "w", newline="") as fh:
    wr = csv.DictWriter(fh, fieldnames=list(rows[1].keys()))
    wr.writeheader()
    wr.writerows(rows)
write(HERE / "anim.extxyz", images)

col = {k: np.array([r.get(k, np.nan) for r in rows], float) for k in rows[1]}
bound = (t_all >= 3.5) & (t_all < 6.0)
qc = {
    "frames": NFRAMES, "fps": FPS, "duration_s": T, "active_li": LI, "N-N_A": round(D_NN, 4),
    "min_N_framework_A": float(np.nanmin(col["min_N_framework_A"])),
    "min_N_other_Li_A": float(np.nanmin(col["min_N_other_Li_A"])),
    "min_Li-N_A": float(np.nanmin(col["Li-N_A"])),
    "j_near_other_Li": int(j_near), "bound_Li-N_range_A": [float(col["Li-N_A"][bound].min()), float(col["Li-N_A"][bound].max())],
    "bound_Li-N-N_range_deg": [float(col["Li-N-N_deg"][bound].min()), float(col["Li-N-N_deg"][bound].max())],
    "glow_frames": [int(i) for i in np.where(col["glow"] == 1)[0][[0, -1]]],
    "glow_seconds": [round(float(t_all[col["glow"] == 1][0]), 3), round(float(t_all[col["glow"] == 1][-1]), 3)],
    "max_step_N_A_per_frame": float(np.nanmax(col["step_N_max_A"])),
    "max_rot_deg_per_frame": float(np.nanmax(col["step_rot_deg"])),
    "loop_seam_A": round(seam, 6), "tumble_search": TUMBLE_SEARCH,
    "limits": {"N-framework": MIN_N_FRAMEWORK, "N-other Li": MIN_N_OTHER_LI, "N-active Li": MIN_N_ACTIVE_LI},
}
(HERE / "qc_summary.json").write_text(json.dumps(qc, indent=1))
print(json.dumps(qc, indent=1))
bad = [r["frame"] for r in rows if r["min_N_framework_A"] < MIN_N_FRAMEWORK
       or r["min_N_other_Li_A"] < MIN_N_OTHER_LI or r["Li-N_A"] < MIN_N_ACTIVE_LI]
assert not bad, f"distance limits broken in frames {bad}"
print("all", NFRAMES, "frames pass the distance checks")
