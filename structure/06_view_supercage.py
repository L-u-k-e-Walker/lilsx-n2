"""06_view_supercage.py - OVITO check renders of the Li-LSX supercage (env viz).

Renders (structure/views/):
  cell_overview.png     whole 672-atom cell, framework sticks + Li spheres
  supercage_cutaway.png one supercage (centre 3/8,3/8,3/8), camera on the [111]
                        12-ring window axis; atoms within 9.5 A of the centre and
                        behind a plane 1.5 A in front of the centre (front wall cut away)
Li colours: SI' grey (hidden in sodalite cages), SII purple, SIII orange.
"""
import warnings
warnings.filterwarnings("ignore", message=".*OVITO.*PyPI")

from pathlib import Path

import numpy as np
from ovito.io import import_file
from ovito.modifiers import CreateBondsModifier, ExpressionSelectionModifier, DeleteSelectedModifier
from ovito.vis import Viewport, TachyonRenderer, BondsVis

HERE = Path(__file__).parent
OUT = HERE / "views"
OUT.mkdir(exist_ok=True)
A = 24.7638
CENTRE = np.array([3/8, 3/8, 3/8]) * A
SITE_COLOURS = {"SI'": (0.55, 0.55, 0.55), "SII": (0.55, 0.25, 0.80), "SIII": (0.95, 0.55, 0.15)}
# atom order in the CIF: framework 0-575, Li SI' 576-607, SII 608-639, SIII 640-671
SITE_OF = {**{i: "SI'" for i in range(576, 608)}, **{i: "SII" for i in range(608, 640)},
           **{i: "SIII" for i in range(640, 672)}}


def style(frame, data):
    types = data.particles.particle_types
    names = {t.id: t.name for t in types.types}
    col = data.particles_.create_property("Color")
    rad = data.particles_.create_property("Radius")
    for i, t in enumerate(types[...]):
        n = names[t]
        if n == "Li":
            col[i] = SITE_COLOURS[SITE_OF[i]]
            rad[i] = 0.75
        elif n == "O":
            col[i] = (0.80, 0.35, 0.30)
            rad[i] = 0.18
        elif n == "Si":
            col[i] = (0.75, 0.70, 0.55)
            rad[i] = 0.22
        else:  # Al
            col[i] = (0.60, 0.70, 0.80)
            rad[i] = 0.22


def pipeline(radius=None, view=None, front_cut=None):
    p = import_file(str(HERE / "LiLSX_P1_MACE.cif"))
    p.modifiers.append(style)
    if radius is not None:
        # keep atoms within `radius` of the supercage centre (minimum image)
        def cut(frame, data):
            pos = data.particles.positions[...]
            d = pos - CENTRE
            d -= np.round(d / A) * A
            sel = data.particles_.create_property("Selection")
            drop = np.linalg.norm(d, axis=1) > radius
            if view is not None:
                # remove the front wall: atoms nearer the camera than front_cut
                drop |= (d @ view) < -front_cut
            sel[...] = drop.astype(int)
            # unwrap kept atoms around the centre so the cage is not split by the cell boundary
            data.particles_.positions_[...] = CENTRE + d
        p.modifiers.append(cut)
        p.modifiers.append(DeleteSelectedModifier())
    cb = CreateBondsModifier(mode=CreateBondsModifier.Mode.Pairwise)
    for t in ("Si", "Al"):
        cb.set_pairwise_cutoff(t, "O", 2.0)
    cb.vis.width = 0.22
    cb.vis.coloring_mode = BondsVis.ColoringMode.ByParticle
    p.modifiers.append(cb)
    p.compute()
    return p


renderer = TachyonRenderer(ambient_occlusion=True, ambient_occlusion_brightness=0.8, antialiasing_samples=12)

# overview
p = pipeline()
p.add_to_scene()
p.source.data.cell.vis.enabled = True
vp = Viewport(type=Viewport.Type.Perspective, camera_dir=(-1, -0.75, -0.55), fov=np.radians(35))
vp.zoom_all()
vp.render_image(size=(960, 720), filename=str(OUT / "cell_overview.png"), background=(1, 1, 1), renderer=renderer)
p.remove_from_scene()

# one supercage, camera on the [111] window axis, front wall cut away
d = np.array([1.0, 1.0, 1.0]) / np.sqrt(3)          # viewing direction (camera -> centre)
p = pipeline(radius=9.5, view=d, front_cut=1.5)
p.add_to_scene()
p.source.data.cell.vis.enabled = False
vp = Viewport(type=Viewport.Type.Perspective, camera_pos=tuple(CENTRE - 45 * d), camera_dir=tuple(d),
              camera_up=(0, 0, 1), fov=np.radians(28))
vp.render_image(size=(960, 720), filename=str(OUT / "supercage_cutaway.png"), background=(1, 1, 1), renderer=renderer)
print("wrote", sorted(x.name for x in OUT.iterdir()))
