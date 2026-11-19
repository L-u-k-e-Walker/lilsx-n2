"""preview.py - quick DRAFT preview of anim.extxyz (env viz), for checking
motion only: 480x360, Tachyon without ambient occlusion, no Li...N glow,
Phase-2 cutaway camera. Frames -> frames/ (gitignored), GIF -> preview_draft.gif.
The real render is Phase 4.
"""
import warnings
warnings.filterwarnings("ignore", message=".*OVITO.*PyPI")
import json
import subprocess
from pathlib import Path

import numpy as np
from ovito.io import import_file
from ovito.modifiers import CreateBondsModifier
from ovito.vis import Viewport, TachyonRenderer, BondsVis

HERE = Path(__file__).parent
FR = HERE / "frames"
FR.mkdir(exist_ok=True)
for f in FR.glob("*.png"):
    f.unlink()
qc = json.loads((HERE / "qc_summary.json").read_text())
LI = qc["active_li"]

p = import_file(str(HERE / "anim.extxyz"))
data0 = p.source.data
A = data0.cell[0, 0]
CENTRE = np.full(3, 3 / 8) * A
VIEW = np.ones(3) / np.sqrt(3)
COL = {"O": (0.80, 0.35, 0.30), "Si": (0.75, 0.70, 0.55), "Al": (0.60, 0.70, 0.80),
       "Li": (0.55, 0.25, 0.80), "N": (0.15, 0.45, 1.00)}
RAD = {"O": 0.18, "Si": 0.22, "Al": 0.22, "Li": 0.70, "N": 0.55}


def cut_and_style(frame, data):
    data.cell_.vis.enabled = False
    pos = data.particles.positions[...]
    d = pos - CENTRE
    d -= np.round(d / A) * A
    names = {t.id: t.name for t in data.particles.particle_types.types}
    el = np.array([names[t] for t in data.particles.particle_types[...]])
    drop = ((np.linalg.norm(d, axis=1) > 9.5) | ((d @ VIEW) < -1.5)) & (el != "N")
    data.particles_.positions_[...] = CENTRE + d
    col = data.particles_.create_property("Color")
    rad = data.particles_.create_property("Radius")
    for e in COL:
        m = el == e
        col[m] = COL[e]
        rad[m] = RAD[e]
    col[LI] = (0.80, 0.45, 1.00)
    rad[LI] = 0.85
    data.particles_.delete_elements(drop)


p.modifiers.append(cut_and_style)
cb = CreateBondsModifier(mode=CreateBondsModifier.Mode.Pairwise)
for t in ("Si", "Al"):
    cb.set_pairwise_cutoff(t, "O", 2.0)
cb.set_pairwise_cutoff("N", "N", 1.3)
cb.vis.width = 0.22
cb.vis.coloring_mode = BondsVis.ColoringMode.ByParticle
p.modifiers.append(cb)
p.add_to_scene()
vp = Viewport(type=Viewport.Type.Perspective, camera_pos=tuple(CENTRE - 45 * VIEW), camera_dir=tuple(VIEW),
              camera_up=(0, 0, 1), fov=np.radians(28))
vp.render_anim(filename=str(FR / "draft_.png"), size=(480, 360), background=(1, 1, 1),
               renderer=TachyonRenderer(ambient_occlusion=False, antialiasing_samples=4), fps=15)
pngs = sorted(FR.glob("draft_*.png"))
print(len(pngs), "frames")
pal = HERE / "frames" / "palette.png"
subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-framerate", "15", "-i", str(FR / "draft_%04d.png"),
                "-vf", "palettegen=stats_mode=diff", str(pal)], check=True)
subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-framerate", "15", "-i", str(FR / "draft_%04d.png"),
                "-i", str(pal), "-lavfi", "paletteuse=dither=sierra2_4a", "-loop", "0",
                str(HERE / "preview_draft.gif")], check=True)
print("wrote preview_draft.gif", (HERE / "preview_draft.gif").stat().st_size // 1024, "KB")
