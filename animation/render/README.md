# Phase 4: render (OVITO 3.17.2 Python module, Tachyon, env `viz`)

Outputs:
* `n2_lilsx.gif`: 960×720, 120 frames, 15 fps, 8.0 s, loops, 0.63 MB (ffmpeg palettegen/paletteuse)
* `n2_lilsx.mp4`: H.264, yuv420p, CRF 18, 0.11 MB (plays more smoothly in PowerPoint than the GIF)
* `still_070.png`: docked frame; `still_020.png`: free-flight frame
* `render_meta.json`: camera, atom counts, timings and file sizes

## Re-render

```bash
conda run -n viz python -I render.py --still 70        # one frame -> still_070.png (~1 s)
conda run -n viz python -I render.py --all             # frames/ + GIF + MP4 (~80 s on 20 cores)
# options: --bg white  --orbit DEG  --size WxH  --out name.png (stills)
```

The trajectory comes from `../trajectory/anim.extxyz`; regenerate it with
`../trajectory/make_trajectory.py`. Frames go to `frames/` (gitignored).
`viz` needs `LD_LIBRARY_PATH=$CONDA_PREFIX/lib`, which is already set as an
env config var. Without it, the pip OVITO's headless Qt plugin can't find
`libOpenGL.so.0`; conda-forge `libopengl`, `libegl` and `libgl` are installed
in the env for this.

## Scene choices

* **Camera.** Fixed, perspective, 22.6° field of view.
  * It starts from the supercage's [111] 12-ring window axis, made
    perpendicular to the docked N₂ axis so Li···N≡N reads side-on.
  * It is then orbited 25° about the N₂ axis for an oblique view into the cage
    (`--orbit`). −30° and 25° were the best of −30/0/25/40/55°; 55° looks
    straight down a cube axis.
  * The aim point is 30 % of the way from the cage centre to the docked N₂.
* **Cutaway.** Framework more than 1.5 Å in front of the cage centre (along
  the view) is not drawn.
  * Whole TO₄ tetrahedra are drawn for T atoms within 13 Å of the cage
    centre, taken from a 3×3×3 block of periodic images.
  * A Li is drawn only if every O within 2.5 Å of it is drawn, so no cation
    floats in space.
  * A depth cue blends far atoms toward the background (up to 62 %).
  * Per frame this draws 32 Si, 34 Al, 165 O and 19 Li, plus the N₂.
* **Styling.**
  * Framework: thin muted sticks (radius 0.09 Å, bond width 0.18 Å).
  * Li⁺: purple spheres (0.55 Å); the binding Li (SIII 662) is larger
    (0.68 Å) and pale violet.
  * N₂: green (0.46 Å).
  * Background: dark.
* **Li···N connector.**
  * It is a viewport overlay (Qt painter) drawn across the gap between the Li
    and N sphere surfaces: one thin (1.4 px), semi-transparent, light-grey
    dashed stroke, with no glow or halo.
  * It shows only while Li–N < 2.8 Å, fading in from 2.8 to 2.12 Å (about
    4 frames each way). Bound Li–N is 2.06–2.10 Å, so it is at full strength
    while docked.
  * No Li–N bond exists in the scene: bonds are made only for Si–O, Al–O and
    N–N.
* **Legend.** An overlay, bottom right: Si, Al, O, Li⁺ (extra-framework), Li⁺
  binding site (SIII), N₂, and "Li⁺···N 2.06 Å (equilibrium)".
  * That distance is the docked MACE-MP-0 + D3 equilibrium Li–N, read from
    `../binding/binding_summary.json`.
  * The 2.8 Å display threshold is internal and not shown.
* **Renderer.** Tachyon with ambient occlusion (12 samples), 16× antialiasing.
  About 0.7 s per frame.
* **GIF.** `palettegen=stats_mode=full`, then
  `paletteuse=dither=sierra2_4a:diff_mode=rectangle`.
