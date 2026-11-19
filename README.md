# lilsx-n2

N₂ adsorption in Li-LSX (Li-exchanged low-silica X, FAU, Si/Al = 1).

| Folder | Status |
|---|---|
| `structure/` | Li-LSX unit cell: built from IZA FAU + literature Li sites and relaxed with MACE-MP-0 + D3, because no open experimental CIF exists. Outputs `structure/LiLSX_P1_MACE.cif` (simulation cell) and `structure/LiLSX_Fd-3_averaged.cif`; method, sources and validation in `structure/README.md`. Done 2026-11-18 |
| `animation/` | 8 s slide GIF/MP4 of N₂ docking end-on at an SIII Li⁺ (scripted illustration; docked geometry from MACE-MP-0 + D3: Li–N 2.063 Å, Li–N–N 171°). Three versions (dark, white, white with an out-of-focus wider framework) in `animation/render/`; pipeline and re-render steps in `animation/README.md`. Done 2026-11-18 |

Environments: `viz` (OVITO 3.17.2, ase, ffmpeg) for rendering; `ml` (ASE 3.29, pymatgen, mace-torch 0.3.17) for the structure work.
Bulky outputs (`*.traj`, `*.out`) are gitignored.
