"""render.py - render the N2-on-Li-LSX storyboard with OVITO (env viz).

usage:
  python render.py --still 70            one frame -> still_070.png
  python render.py --all                 120 frames -> frames/, then GIF + MP4
  options: --bg dark|white  --size 960x720

Scene
  * camera fixed, looking into the supercage centred at (3/8,3/8,3/8); tilted
    from the [111] window axis so the docked Li...N-N is seen side-on
  * cutaway: framework atoms nearer the camera than a plane CUT A in front of
    the cage centre, or farther than R_KEEP from it, are not drawn
  * framework: thin muted sticks (Si, Al, O); Li+ purple spheres, the active
    Li (SIII 662) larger and brighter; N2 bright blue
  * Li...N interaction: dashed glowing connector painted by a viewport
    overlay across the gap between the Li and N sphere surfaces, only when
    Li-N < 2.8 A (fades in between 2.8 and 2.12 A, ~4 frames). Never a bond.
  * Tachyon ray tracer with ambient occlusion
"""
import warnings
warnings.filterwarnings("ignore", message=".*OVITO.*PyPI")

import argparse
import json
import subprocess
import time
from pathlib import Path

import numpy as np
from ase.io import read
from ovito.io.ase import ase_to_ovito
from ovito.modifiers import CreateBondsModifier
from ovito.pipeline import Pipeline, StaticSource
from ovito.vis import BondsVis, PythonViewportOverlay, TachyonRenderer, Viewport, ViewportOverlayInterface
from PySide6.QtCore import QPointF, Qt
from PySide6.QtGui import QColor, QPen, QRadialGradient, QBrush

HERE = Path(__file__).parent
ROOT = HERE.parents[1]
ap = argparse.ArgumentParser()
ap.add_argument("--still", type=int)
ap.add_argument("--all", action="store_true")
ap.add_argument("--bg", default="dark", choices=["dark", "white"])
ap.add_argument("--size", default="960x720")
ap.add_argument("--orbit", type=float, default=None, help="camera orbit about the bound N2 axis (deg)")
ap.add_argument("--out", default=None, help="still output file name")
args = ap.parse_args()
W, H = map(int, args.size.split("x"))

# ---------------------------------------------------------------- look
GLOW_ON, GLOW_FULL = 2.8, 2.12        # A: connector appears below 2.8, full by 2.12 (bound: 2.06-2.10)
CUT, R_T = 1.5, 13.0                  # A: cutaway plane offset; T atoms kept within R_T of the cage centre
ORBIT_DEG = 25.0 if args.orbit is None else args.orbit   # camera rotation about the bound N2 axis
FADE_START, FADE_LEN, FADE_MAX = 2.0, 14.0, 0.62          # depth cue: blend far atoms towards the background
if args.bg == "dark":
    BG = (0.035, 0.045, 0.075)
    FW = {"Si": (0.62, 0.60, 0.55), "Al": (0.50, 0.58, 0.68), "O": (0.66, 0.42, 0.40)}
    GLOW_RGB = (255, 214, 102)
else:
    BG = (1.0, 1.0, 1.0)
    FW = {"Si": (0.70, 0.66, 0.56), "Al": (0.55, 0.63, 0.74), "O": (0.78, 0.46, 0.42)}
    GLOW_RGB = (235, 150, 20)
LI_RGB, LI_ACTIVE_RGB = (0.60, 0.38, 0.90), (0.86, 0.62, 1.00)
N_RGB = (0.22, 0.60, 1.00)
R_FW, BOND_W = 0.09, 0.18
R_LI, R_LI_ACTIVE, R_N = 0.55, 0.68, 0.46

# ---------------------------------------------------------------- data
traj = HERE.parent / "trajectory"
frames = read(traj / "anim.extxyz", index=":")
qc = json.loads((traj / "qc_summary.json").read_text())
LI = qc["active_li"]
NFW = len(frames[0]) - 2
A = frames[0].cell[0, 0]
CENTRE = np.full(3, 3 / 8) * A
VIEW0 = np.ones(3) / np.sqrt(3)


def mic(v):
    return v - np.round(v / A) * A


# bound N2 axis (middle of the docked phase) -> camera direction
mid = frames[70]
nin = mid.info["n_inner"]
nout = NFW + 1 if nin == NFW else NFW
u_b = mic(mid.positions[nout] - mid.positions[nin])
u_b /= np.linalg.norm(u_b)


def rotate(v, axis, ang):
    axis = axis / np.linalg.norm(axis)
    return (v * np.cos(ang) + np.cross(axis, v) * np.sin(ang) + axis * (axis @ v) * (1 - np.cos(ang)))


# start from the [111] window axis made perpendicular to the N2 axis (N2 side-on),
# then orbit about the N2 axis: stays side-on, gives an oblique view into the cage
DIR = VIEW0 - (VIEW0 @ u_b) * u_b
DIR /= np.linalg.norm(DIR)
DIR = rotate(DIR, u_b, np.radians(ORBIT_DEG))
tilt_deg = float(np.degrees(np.arccos(DIR @ VIEW0)))
com_b = CENTRE + mic(0.5 * (mid.positions[nin] + mid.positions[nout]) - CENTRE)
CAM_DIST = 60.0
_up = np.array([0, 0, 1.0]) - DIR[2] * DIR
_up /= np.linalg.norm(_up)
_right = np.cross(DIR, _up)
# aim between the cage centre and the docked N2; vertical half-height 0.20 (by eye)
TARGET = 0.7 * CENTRE + 0.3 * com_b
FOV = 2 * np.arctan(0.20)

# framework selection from a 3x3x3 block of periodic images (R_T > a/2):
# whole TO4 tetrahedra for T atoms inside the region; Li only when every O it
# touches (< 2.5 A) is drawn, so no Li floats in empty space
fw0 = frames[0][:NFW]
sym0 = np.array(fw0.get_chemical_symbols())
shifts = np.array([[i, j, k] for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1)]) * A
sup_pos = (fw0.positions[None, :, :] + shifts[:, None, :]).reshape(-1, 3)
sup_idx = np.tile(np.arange(NFW), len(shifts))
sup_sym = sym0[sup_idx]
dd = sup_pos - CENTRE
inside = (np.linalg.norm(dd, axis=1) <= R_T + 3.0) & ((dd @ DIR) >= -CUT - 3.0)
cand = np.where(inside)[0]
cpos, csym = sup_pos[cand], sup_sym[cand]
cd = cpos - CENTRE
isT = np.isin(csym, ["Si", "Al"])
isO, isLi = csym == "O", csym == "Li"
keepT = isT & (np.linalg.norm(cd, axis=1) <= R_T) & ((cd @ DIR) >= -CUT)
DTO = np.linalg.norm(cpos[:, None, :] - cpos[None, keepT, :], axis=2)
keepO = isO & (DTO < 2.0).any(axis=1)
keep_c = keepT | keepO
DLiO = np.linalg.norm(cpos[isLi][:, None, :] - cpos[isO][None, :, :], axis=2)
li_rows = np.where(isLi)[0]
o_rows = np.where(isO)[0]
keepLi = np.zeros(len(cand), bool)
for r, row in enumerate(li_rows):
    if (cd[row] @ DIR) < -CUT or np.linalg.norm(cd[row]) > R_T + 1.0:
        continue
    neigh = o_rows[DLiO[r] < 2.5]
    keepLi[row] = len(neigh) > 0 and keep_c[neigh].all()
# the active Li: its image nearest the cage centre
act_rows = np.where(sup_idx[cand] == LI)[0]
act_row = act_rows[np.argmin(np.linalg.norm(cd[act_rows], axis=1))]
keepLi[act_row] = True
keep_rows = np.where(keep_c | keepLi)[0]
fw_pos = cpos[keep_rows]
fw_sym = csym[keep_rows]
active_k = int(np.where(keep_rows == act_row)[0][0])
depth = (fw_pos - TARGET) @ DIR
fade = np.clip((depth - FADE_START) / FADE_LEN, 0, FADE_MAX)


def build(fr):
    from ase import Atoms
    at = Atoms(list(fw_sym), positions=fw_pos)
    n = fr[NFW:NFW + 2].copy()
    n.positions = CENTRE + mic(fr.positions[NFW:NFW + 2] - CENTRE)
    at.extend(n)
    at.pbc = False
    data = ase_to_ovito(at)
    data.cell_.vis.enabled = False
    col = data.particles_.create_property("Color")
    rad = data.particles_.create_property("Radius")
    sym = np.array(at.get_chemical_symbols())
    for e, c in FW.items():
        col[sym == e] = c
        rad[sym == e] = R_FW
    col[sym == "Li"] = LI_RGB
    rad[sym == "Li"] = R_LI
    col[active_k] = LI_ACTIVE_RGB
    rad[active_k] = R_LI_ACTIVE
    col[sym == "N"] = N_RGB
    rad[sym == "N"] = R_N
    # depth cue on the static framework + spectator Li (not the active Li, not N2)
    nf = len(fw_sym)
    f = fade.copy()
    f[active_k] = 0.0
    c = np.array(col[:nf])
    col[:nf] = (1 - f[:, None]) * c + f[:, None] * np.array(BG)
    return data, at


# ---------------------------------------------------------------- glow overlay
class Glow(ViewportOverlayInterface):
    """dashed glowing Li...N connector between the sphere surfaces"""
    state = {"on": 0.0}

    def render(self, canvas, **kwargs):
        s = self.state
        if s["on"] <= 0:
            return
        p_li = canvas.project_location(tuple(s["li"]))
        p_n = canvas.project_location(tuple(s["n"]))
        if p_li is None or p_n is None:
            return
        Wc, Hc = canvas.logical_size

        def px(p):
            return np.array([p[0] * Wc, (1 - p[1]) * Hc])

        a, b = px(p_li), px(p_n)
        v = b - a
        L = np.linalg.norm(v)
        if L < 1:
            return
        v /= L
        r_li = canvas.project_length(tuple(s["li"]), R_LI_ACTIVE) * Hc
        r_n = canvas.project_length(tuple(s["n"]), R_N) * Hc
        a2, b2 = a + v * r_li * 0.92, b - v * r_n * 0.92
        if np.linalg.norm(b2 - a2) < 2 or (b2 - a2) @ v <= 0:
            return
        f = s["on"]
        R, G, B = GLOW_RGB
        with canvas.qt_painter() as qp:
            qp.setRenderHint(qp.RenderHint.Antialiasing)
            # soft halo around the gap
            c = 0.5 * (a2 + b2)
            rad = 0.9 * np.linalg.norm(b2 - a2) + 10
            g = QRadialGradient(QPointF(*c), rad)
            g.setColorAt(0.0, QColor(R, G, B, int(110 * f)))
            g.setColorAt(1.0, QColor(R, G, B, 0))
            qp.setPen(Qt.NoPen)
            qp.setBrush(QBrush(g))
            qp.drawEllipse(QPointF(*c), rad, rad)
            # dashed line with layered glow strokes (wide/faint -> thin/bright)
            seg = np.linalg.norm(b2 - a2)
            for width, alpha in ((9, 40), (6, 70), (3.5, 140), (1.8, 255)):
                pen = QPen(QColor(R, G, B, int(alpha * f)))
                pen.setWidthF(width)
                pen.setCapStyle(Qt.RoundCap)
                # dash pattern in units of pen width; ~3 dashes across the gap
                dash = max(seg / 5.5, 4.0)
                pen.setDashPattern([dash / width, dash * 0.8 / width])
                qp.setPen(pen)
                qp.drawLine(QPointF(*a2), QPointF(*b2))


glow = Glow()


class Legend(ViewportOverlayInterface):
    """atom key, bottom-right: swatches painted with Qt, text by OVITO"""
    ITEMS = [("sphere", FW["Si"], 0.45, "Si"), ("sphere", FW["Al"], 0.45, "Al"), ("sphere", FW["O"], 0.45, "O"),
             ("sphere", LI_RGB, 0.75, "Li⁺ (extra-framework)"),
             ("sphere", LI_ACTIVE_RGB, 0.9, "Li⁺ binding site (SIII)"),
             ("n2", N_RGB, 0.7, "N₂"),
             ("dash", None, 0.0, "Li⁺···N interaction (< 2.8 Å)")]

    def render(self, canvas, **kwargs):
        Wc, Hc = canvas.logical_size
        fs = 0.026                                      # font size (fraction of height)
        row = 0.040 * Hc                                # row pitch, px
        box_w, box_h = 0.335 * Wc, row * (len(self.ITEMS) + 0.6)
        x0, y0 = Wc - box_w - 0.025 * Wc, Hc - box_h - 0.03 * Hc     # top-left of box (px, Qt coords)
        text_rgb = (0.86, 0.88, 0.92) if args.bg == "dark" else (0.15, 0.17, 0.22)
        with canvas.qt_painter() as qp:
            qp.setRenderHint(qp.RenderHint.Antialiasing)
            bgc = QColor(int(BG[0] * 255), int(BG[1] * 255), int(BG[2] * 255), 205)
            edge = QColor(255, 255, 255, 40) if args.bg == "dark" else QColor(0, 0, 0, 40)
            qp.setBrush(QBrush(bgc))
            qp.setPen(QPen(edge, 1.0))
            qp.drawRoundedRect(x0, y0, box_w, box_h, 8, 8)
            for k, (kind, rgb, size, _) in enumerate(self.ITEMS):
                cx, cy = x0 + 0.075 * box_w, y0 + row * (k + 0.8)
                r = 0.42 * row * size
                if kind == "sphere":
                    g = QRadialGradient(QPointF(cx - 0.35 * r, cy - 0.35 * r), 1.4 * r)
                    c = QColor(int(rgb[0] * 255), int(rgb[1] * 255), int(rgb[2] * 255))
                    g.setColorAt(0.0, c.lighter(135))
                    g.setColorAt(1.0, c.darker(150))
                    qp.setPen(Qt.NoPen)
                    qp.setBrush(QBrush(g))
                    qp.drawEllipse(QPointF(cx, cy), r, r)
                elif kind == "n2":
                    c = QColor(int(rgb[0] * 255), int(rgb[1] * 255), int(rgb[2] * 255))
                    qp.setPen(Qt.NoPen)
                    qp.setBrush(QBrush(c))
                    qp.drawEllipse(QPointF(cx - 0.55 * r, cy), 0.8 * r, 0.8 * r)
                    qp.drawEllipse(QPointF(cx + 0.55 * r, cy), 0.8 * r, 0.8 * r)
                else:
                    R, G, B = GLOW_RGB
                    for width, alpha in ((6, 60), (2.2, 255)):
                        pen = QPen(QColor(R, G, B, alpha))
                        pen.setWidthF(width)
                        pen.setCapStyle(Qt.RoundCap)
                        pen.setDashPattern([3.0 / width * 2, 2.6 / width * 2])
                        qp.setPen(pen)
                        qp.drawLine(QPointF(cx - 0.36 * row, cy), QPointF(cx + 0.36 * row, cy))
        for k, (_, _, _, label) in enumerate(self.ITEMS):
            ty = y0 + row * (k + 0.8)
            canvas.draw_text(label, pos=((x0 + 0.15 * box_w) / Wc, 1 - ty / Hc), font_size=fs,
                             anchor="west", color=text_rgb)


legend = Legend()

# ---------------------------------------------------------------- scene
data0, _ = build(frames[0])
pipe = Pipeline(source=StaticSource(data=data0))
cb = CreateBondsModifier(mode=CreateBondsModifier.Mode.Pairwise)
for t in ("Si", "Al"):
    cb.set_pairwise_cutoff(t, "O", 2.0)
cb.set_pairwise_cutoff("N", "N", 1.3)          # N2 only; no Li-N bond is ever created
cb.vis.width = BOND_W
cb.vis.coloring_mode = BondsVis.ColoringMode.ByParticle
pipe.modifiers.append(cb)
pipe.add_to_scene()
vp = Viewport(type=Viewport.Type.Perspective, camera_pos=tuple(TARGET - CAM_DIST * DIR),
              camera_dir=tuple(DIR), camera_up=(0, 0, 1), fov=FOV)
vp.overlays.append(PythonViewportOverlay(delegate=glow))
vp.overlays.append(PythonViewportOverlay(delegate=legend))
renderer = TachyonRenderer(ambient_occlusion=True, ambient_occlusion_brightness=0.8,
                           ambient_occlusion_samples=12, antialiasing_samples=16,
                           direct_light_intensity=0.95)


def glow_factor(dist):
    x = np.clip((GLOW_ON - dist) / (GLOW_ON - GLOW_FULL), 0, 1)
    return float(x * x * (3 - 2 * x)) if dist < GLOW_ON else 0.0


def render_frame(i, path):
    fr = frames[i]
    data, at = build(fr)
    pipe.source.data = data
    n_k = len(at) - 2 + (fr.info["n_inner"] - NFW)
    glow.state = {"on": glow_factor(fr.info["li_n_A"]), "li": at.positions[active_k], "n": at.positions[n_k]}
    vp.render_image(size=(W, H), filename=str(path), background=BG, renderer=renderer)
    return glow.state["on"]


meta = {"camera_dir": DIR.round(4).tolist(), "tilt_from_111_deg": round(tilt_deg, 1),
        "target": TARGET.round(3).tolist(), "camera_distance_A": CAM_DIST, "fov_deg": round(np.degrees(FOV), 2),
        "orbit_deg": ORBIT_DEG, "atoms_drawn": {e: int((fw_sym == e).sum()) for e in ("Si", "Al", "O", "Li")}, "bg": args.bg, "size": [W, H]}
print(json.dumps(meta))

if args.still is not None:
    t0 = time.time()
    out = HERE / (args.out or f"still_{args.still:03d}{'' if args.bg == 'dark' else '_white'}.png")
    f = render_frame(args.still, out)
    print(f"wrote {out.name}  glow={f:.2f}  Li-N={frames[args.still].info['li_n_A']:.3f} A  "
          f"{time.time() - t0:.1f} s")

if args.all:
    FR = HERE / "frames"
    FR.mkdir(exist_ok=True)
    for p in FR.glob("*.png"):
        p.unlink()
    t0 = time.time()
    for i in range(len(frames)):
        render_frame(i, FR / f"f_{i:04d}.png")
    t_render = time.time() - t0
    suffix = "" if args.bg == "dark" else "_white"
    gif, mp4, pal = HERE / f"n2_lilsx{suffix}.gif", HERE / f"n2_lilsx{suffix}.mp4", FR / "palette.png"
    subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-framerate", "15", "-i", str(FR / "f_%04d.png"),
                    "-vf", "palettegen=max_colors=256:stats_mode=full", str(pal)], check=True)
    subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-framerate", "15", "-i", str(FR / "f_%04d.png"),
                    "-i", str(pal), "-lavfi", "paletteuse=dither=sierra2_4a:diff_mode=rectangle",
                    "-loop", "0", str(gif)], check=True)
    subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-framerate", "15", "-i", str(FR / "f_%04d.png"),
                    "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "18", "-preset", "slow",
                    "-movflags", "+faststart", str(mp4)], check=True)
    meta.update({"render_s": round(t_render, 1), "gif_MB": round(gif.stat().st_size / 1e6, 2),
                 "mp4_MB": round(mp4.stat().st_size / 1e6, 2)})
    (HERE / f"render_meta{suffix}.json").write_text(json.dumps(meta, indent=1))
    print(json.dumps(meta, indent=1))
