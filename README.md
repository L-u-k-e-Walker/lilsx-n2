# lilsx-n2

N₂ adsorption in Li-LSX (Li-exchanged low-silica X, FAU, Si/Al = 1).

| Folder | Status |
|---|---|
| `structure/` | Li-LSX unit cell: built from IZA FAU + literature Li sites and relaxed with MACE-MP-0 + D3, because no open experimental CIF exists. Outputs `structure/LiLSX_P1_MACE.cif` (simulation cell) and `structure/LiLSX_Fd-3_averaged.cif`; method, sources and validation in `structure/README.md`. Done 2026-11-18 |
| `animation/binding/` | Phase 2 of the N₂ adsorption animation: N₂ end-on at SIII Li 662, MACE-MP-0 + D3 (N₂ only relaxed): Li–N 2.063 Å, Li–N–N 171°; compared with gas-phase Li⁺–N₂ ab initio (≈ 2.07 Å). See `animation/binding/README.md` |

Environments: `viz` (OVITO 3.17.2, ase, ffmpeg) for rendering; `ml` (ASE 3.29, pymatgen, mace-torch 0.3.17) for the structure work.
Bulky outputs (`*.traj`, `*.out`) are gitignored.
