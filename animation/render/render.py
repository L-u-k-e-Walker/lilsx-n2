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
    Li-N < 2.8 A (fades in between 2.8 and 2.45 A). Never a bond.
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
args = ap.parse_args()
W, H = map(int, args.size.split("x"))

# ---------------------------------------------------------------- look
GLOW_ON, GLOW_FULL = 2.8, 2.45        # A: connector appears below 2.8, full by 2.45
CUT, R_KEEP = 1.5, 10.0               # A: cutaway plane offset, keep radius
TILT_FRACTION = 1.0                   # 1 = camera exactly perpendicular to the bound N2 axis
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
li_pos = CENTRE + mic(mid.positions[LI] - CENTRE)
nin = mid.info["n_inner"]
nout = NFW + 1 if nin == NFW else NFW
u_b = mic(mid.positions[nout] - mid.positions[nin])
u_b /= np.linalg.norm(u_b)
DIR = VIEW0 - TILT_FRACTION * (VIEW0 @ u_b) * u_b
DIR /= np.linalg.norm(DIR)
tilt_deg = float(np.degrees(np.arccos(DIR @ VIEW0)))
com_b = CENTRE + mic(0.5 * (mid.positions[nin] + mid.positions[nout]) - CENTRE)
# framing: aim slightly lower-left of the cage centre (where the N2 path is),
# vertical half-height 0.19 in view-angle units (set by eye): whole 12-ring nearly in frame,
# N2 path ~30 % of the image height
CAM_DIST = 60.0
_up = np.array([0, 0, 1.0]) - DIR[2] * DIR
_up /= np.linalg.norm(_up)
_right = np.cross(DIR, _up)
AIM_OFFSET = (-0.02, 0.0)
TARGET = CENTRE + CAM_DIST * (AIM_OFFSET[0] * _right + AIM_OFFSET[1] * _up)
FOV = 2 * np.arctan(0.19)

# static framework subset (cutaway), unwrapped around the cage centre
fw0 = frames[0][:NFW]
d = mic(fw0.positions - CENTRE)
keep = (np.linalg.norm(d, axis=1) <= R_KEEP) & ((d @ DIR) >= -CUT)
keep_idx = np.where(keep)[0]
fw_pos = CENTRE + d[keep_idx]
fw_sym = np.array(fw0.get_chemical_symbols())[keep_idx]
active_k = int(np.where(keep_idx == LI)[0][0])


def build(fr):
    at = fw0[keep_idx].copy()
    at.positions = fw_pos
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
        "atoms_drawn_framework": int(len(keep_idx)), "bg": args.bg, "size": [W, H]}
print(json.dumps(meta))

if args.still is not None:
    t0 = time.time()
    out = HERE / f"still_{args.still:03d}{'' if args.bg == 'dark' else '_white'}.png"
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
