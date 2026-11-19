# Phase 3: scripted N₂ adsorption trajectory (storyboard, not MD)

`make_trajectory.py` (env `viz`) writes `anim.extxyz`. It has 120 frames at
15 fps (8.0 s), each with 674 atoms: the static Li-LSX cell plus one N₂. The
active Li is SIII Li 662, and the bound pose comes from Phase 2 (Li–N 2.063 Å,
Li–N–N 171°).

| Time | Frames | Story |
|---|---|---|
| 0–2.5 s | 0–37 | N₂ drifts in from near the supercage centre, tumbling (228° of rotation) |
| 2.5–3.5 s | 38–52 | turns end-on and docks; orientation blends from tumble to bound (eased) |
| 3.5–6.0 s | 53–89 | bound: smooth wobble of the Li-side N (±0.06 Å) and libration (±7°) |
| 6.0–8.0 s | 90–119 | lifts off end-on, tilted away from the neighbouring SII Li 617, then drifts back towards the cage centre |

## How the motion is made

* **Path.** The N₂ position follows a C¹ cubic Hermite spline through
  keyframes: start, pre-dock (4.4 Å out from the Li), docked (held 3.5–6.0 s),
  lift-off (3 Å out), and back to the start. Docking and lift-off ease in and
  out because the spline has zero velocity at the docked keyframes.
* **Orientation.** The free tumble is two steady rotations about fixed axes,
  at 2 and 1 turns per 8 s. It blends into the bound orientation (spherical
  interpolation, smootherstep weight) from 1.6 to 3.5 s, and back out from
  6.3 to 7.9 s.
* **Choice of tumble axes.** The axes and starting phase were chosen by a
  search over 300 random sets. The winner gives the slowest worst-case
  rotation step while still tumbling at least 150° in 0–2 s and 90° in
  6.8–8 s, and passing every distance check.
* **Pivot.** The molecule's pivot moves with the same blend weight: the N₂
  centre of mass in free flight, the Li-side N when docked.
* **"Thermal" noise.** Only sums of low-frequency sinusoids, with frequencies
  that are integer multiples of 1/8 Hz and fixed random phases. There is no
  per-frame random jitter.
* **Seamless loop.** Every term has an 8 s period, and the path ends where it
  started with the same velocity. The step from frame 119 to frame 0 is as
  small as any other step (seam error 0.0 Å), so the GIF loops without a
  jump. The trade-off: "drifts away" ends back near the cage centre.

## Checks (`qc_summary.json`, per frame in `trajectory_qc.csv`)

| Check | Limit | Result |
|---|---|---|
| N to framework O/Si/Al, every frame | ≥ 2.9 Å | min 3.21 Å |
| N to any non-active Li | ≥ 2.9 Å | min 2.93 Å (SII Li 617, during lift-off) |
| N to active Li | ≥ 1.95 Å | min 2.06 Å |
| Bound Li–N / Li–N–N | – | 2.06–2.10 Å / 169.5–177.9° |
| Largest N step per frame | smooth | 0.33 Å |
| Largest N₂ rotation per frame | smooth | 8.7° |
| Glow (Li–N < 2.8 Å) | – | frames 46–95 (3.07–6.33 s) |

The script aborts if any frame breaks a distance limit. Each frame's info
line carries `time_s`, `li_n_A`, `li_n_n_deg`, `glow`, `active_li`, `n_inner`
and `w_bind` for the renderer.

## Draft preview

`preview.py` writes `preview_draft.gif` (480×360, no ambient occlusion, no
glow, 72 KB). It is for checking the motion only; the real render is Phase 4.
Frames go to `frames/`, which is gitignored.
