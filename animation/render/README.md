# Render (OVITO 3.17.2 Python module, Tachyon, env `viz`)

Three versions of the same 8 s loop (960×720, 120 frames, 15 fps):

| Variant | GIF | MP4 | Stills | Look |
|---|---|---|---|---|
| `dark` (main) | `n2_lilsx_dark.gif` | `n2_lilsx_dark.mp4` | `still_070_dark.png`, `still_020_dark.png` | dark background, 13 Å framework shell, depth fade |
| `white` | `n2_lilsx_white.gif` | `n2_lilsx_white.mp4` | `still_070_white.png`, `still_020_white.png` | same scene on white |
| `white_dof` | `n2_lilsx_white_dof.gif` | `n2_lilsx_white_dof.mp4` | `still_070_white_dof.png`, `still_020_white_dof.png` | white, wider 19 Å framework shell, depth of field focused on the binding Li |

File sizes and timings are in `render_meta_<variant>.json`. Frame 70 is
docked; frame 20 is free flight. The MP4 (H.264, yuv420p, CRF 18) plays more
smoothly in PowerPoint than the GIF.

## Re-render

```bash
conda run -n viz python -I render.py --still 70 --variant dark     # one frame (~1 s; white_dof ~7 s)
conda run -n viz python -I render.py --all --variant dark          # frames/dark/ + GIF + MP4 (~1.5 min)
conda run -n viz python -I render.py --all --variant white         # ~1.5 min
conda run -n viz python -I render.py --all --variant white_dof     # ~15 min (64 samples/pixel for the blur)
# other options: --orbit DEG (camera angle), --size WxH, --out name.png (stills)
```

The trajectory comes from `../trajectory/anim.extxyz`; regenerate it with
`../trajectory/make_trajectory.py`. Frames go to `frames/<variant>/`, which is
gitignored.

`viz` needs `LD_LIBRARY_PATH=$CONDA_PREFIX/lib`, which is already set as an
env config var. Without it, the pip OVITO's headless Qt plugin can't find
`libOpenGL.so.0`; conda-forge `libopengl`, `libegl` and `libgl` are installed
in the env for this.

## Scene choices

* **Camera.** Fixed, perspective, 22.6° field of view.
  * It starts from the supercage's [111] 12-ring window axis, made
    perpendicular to the docked N₂ axis so Li···N≡N reads side-on.
  * It is then orbited 25° about the N₂ axis for an oblique view into the cage
    (`--orbit`; picked from −30/0/25/40/55°).
  * The aim point is 30 % of the way from the cage centre to the docked N₂.
* **Cutaway.** Framework more than 1.5 Å in front of the cage centre (along
  the view) is not drawn.
  * Whole TO₄ tetrahedra are drawn for T atoms within 13 Å of the cage
    centre (19 Å for `white_dof`), taken from a 3×3×3 block of periodic
    images.
  * A Li is drawn only if every O within 2.5 Å of it is drawn, so no cation
    floats in space.
  * A depth cue blends far atoms toward the background: up to 62 %, or 35 %
    for `white_dof`, where the blur carries the depth.
* **Depth of field (`white_dof`).** Tachyon `depth_of_field`, focal plane at
  the binding Li (camera-to-Li distance along the view), aperture 0.25 Å, 64
  antialiasing samples. The binding site and near cage wall stay sharp; the
  wider framework behind is softly out of focus.
* **Styling.**
  * Framework: muted sticks (radius 0.13 Å, bond width 0.26 Å).
  * Li⁺: purple spheres (0.55 Å); the binding Li (SIII 662) is larger
    (0.68 Å) and pale violet.
  * N₂: green (0.46 Å).
* **Li···N connector.**
  * It is a viewport overlay (Qt painter) drawn across the gap between the Li
    and N sphere surfaces: one thin (1.4 px), semi-transparent dashed stroke
    (light grey on dark, dark grey on white), with no glow or halo.
  * It shows only while Li–N < 2.8 Å, fading in from 2.8 to 2.12 Å (about
    4 frames each way). Bound Li–N is 2.06–2.10 Å.
  * No Li–N bond exists in the scene: bonds are made only for Si–O, Al–O and
    N–N.
* **Legend.** An overlay, bottom right: Si, Al, O, Li⁺ (extra-framework), Li⁺
  binding site (SIII), N₂, and "Li⁺···N 2.06 Å (equilibrium)".
  * That distance is the docked MACE-MP-0 + D3 equilibrium Li–N, read from
    `../binding/binding_summary.json`.
  * The 2.8 Å display threshold is internal and not shown.
* **Renderer.** Tachyon with ambient occlusion (12 samples), 16× antialiasing
  (64× for `white_dof`).
* **GIF.** `palettegen=stats_mode=full`, then
  `paletteuse=dither=sierra2_4a:diff_mode=rectangle`.
