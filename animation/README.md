# N₂ adsorbing at Li⁺ in Li-LSX: 8 s slide animation

One N₂ molecule drifts into a supercage of Li-LSX, docks end-on at an exposed
Li⁺ (SIII site), stays for 2 s, then leaves. The loop is seamless.

**This is a scripted illustration, not a simulation.** Only the docked
geometry comes from a calculation: an N₂-only relaxation with MACE-MP-0 + D3,
framework fixed.

**Final files** (`render/`):

| Version | GIF | MP4 |
|---|---|---|
| dark (main) | `n2_lilsx_dark.gif` | `n2_lilsx_dark.mp4` |
| white | `n2_lilsx_white.gif` | `n2_lilsx_white.mp4` |
| white, wider framework out of focus | `n2_lilsx_white_dof.gif` | `n2_lilsx_white_dof.mp4` |

Each is 960×720, 120 frames, 15 fps, 8.0 s. Use the MP4 in PowerPoint (it
plays more smoothly); the GIF is for anything else.

## Pipeline

| Step | Folder | Env | What |
|---|---|---|---|
| 1 Structure | `../structure/` | `ml` | Li-LSX cell built from IZA FAU (Si/Al = 1, strict alternation) with literature Li sites, relaxed with MACE-MP-0 + D3 (no open experimental CIF exists). `LiLSX_P1_MACE.cif`: a = 24.764 Å, Li₉₆Al₉₆Si₉₆O₃₈₄, charge neutral |
| 2 Binding | `binding/` | `ml` | N₂ end-on at SIII Li 662, N₂ only relaxed: Li–N 2.063 Å, Li–N–N 171°; compared with gas-phase Li⁺–N₂ ab initio (≈ 2.07 Å) |
| 3 Trajectory | `trajectory/` | `viz` | scripted storyboard → `anim.extxyz`: Hermite-spline path, quaternion blending, smooth periodic noise, per-frame overlap checks |
| 4–5 Render | `render/` | `viz` | OVITO + Tachyon (ambient occlusion), cutaway, legend, dashed Li···N connector → GIF (ffmpeg palette method) + H.264 MP4 |

Storyboard: approach and tumble 0–2.5 s; dock 2.5–3.5 s; bound 3.5–5.5 s;
leave 5.5–8.0 s.

## Re-render from scratch

```bash
# 1-2 (only if the structure or binding changes; needs the GPU)
cd ../structure && conda run -n ml python -I 01_build.py && for c in 0 1 2; do conda run -n ml python -I 02_relax.py $c; done && conda run -n ml python -I 03_analyze.py
cd ../animation/binding && conda run -n ml python -I 01_bind_n2.py && conda run -n ml python -I 02_tilt_scan.py && conda run -n viz python -I 03_summary_view.py
# 3 trajectory (~5 s)
cd ../trajectory && conda run -n viz python -I make_trajectory.py
# 4 render (dark/white ~1.5 min each, white_dof ~15 min)
cd ../render && for v in dark white white_dof; do conda run -n viz python -I render.py --all --variant $v; done
```

To change only the look (colours, sizes, camera, legend), edit the
constants at the top of `render/render.py` and re-run step 4. To change only
the timing, edit `T_DOCK`, `T_BOUND_END`, `T_OFF` and `W_OFF` in
`trajectory/make_trajectory.py`, then re-run steps 3 and 4.

## Environment `viz`

Python 3.12 with ase, numpy, ffmpeg and ovito 3.17.2 (pip). The pip OVITO
also needs conda-forge `libopengl`, `libegl` and `libgl`, plus
`LD_LIBRARY_PATH=$CONDA_PREFIX/lib`, which is set as an env config var with
`conda env config vars set`.

Rendered frames (`frames/`) are gitignored everywhere.
