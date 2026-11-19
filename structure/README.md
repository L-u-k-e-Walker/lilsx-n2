# Li-LSX unit cell (FAU, Si/Al = 1, Li⁺ at crystallographic sites)

**Result.** No experimental CIF for dehydrated Li-LSX could be found in open
sources (search 2026-11-18). So the cell was **built** from the IZA FAU framework
with strict Si/Al alternation and 96 Li⁺ at the SI′/SII/SIII sites reported in
the Li-LSX literature, then relaxed with MACE-MP-0 + D3. It was checked against
the two open experimental Li-faujasite neutron structures (Li-X and Li-Y in
COD). Treat it as a **model**, not a measured structure.

| File | What it is | Use it for |
|---|---|---|
| `LiLSX_P1_MACE.cif` | ordered relaxed cell, 672 atoms, P1 listing, cubic a = 24.764 Å | simulations (GCMC, MD, DFT) |
| `LiLSX_Fd-3_averaged.cif` | the same model averaged onto Fd-3 (origin 2): Si, Al, O1–O4, Li SI′ 32e, SII 32e, SIII 48f with occupancy 2/3 | crystallographic description, comparing with refinements |
| `results/sites_Fd-3.csv` | averaged asymmetric unit + rms spread over each orbit | |
| `results/li_sites_P1.csv` | site label (SI′/SII/SIII), coordination and shortest Li–O for every Li in the P1 file | assigning Li types in a force field |
| `results/geometry.json` | all checks below (bond lengths, angles, Li coordination, comparison with experiment, symmetry, CIF round trip) | |
| `results/candidates.csv` | energies of the three SIII arrangements | |
| `results/compare_shepelev1990*.{json,csv}` | site-by-site comparison with the published Si/Al-ordered Li,Na-X structure (section 5) | |

Formula: Li₉₆[Al₉₆Si₉₆O₃₈₄], i.e. one Li per Al. In the P1 CIF the labels are
element symbols only (RASPA reads labels as pseudo-atom types), and the atom
order is framework (0–575), then Li SI′ (576–607), SII (608–639), SIII (640–671).

## 1. Search for an experimental CIF

| Source | Content | Usable? |
|---|---|---|
| Plévert, Di Renzo, Fajula, Chiari, *J. Phys. Chem. B* **101** (1997) 10340, [doi:10.1021/jp9714330](https://pubs.acs.org/doi/10.1021/jp9714330) | Dehydrated Li-LSX, neutron powder diffraction; cubic Fd-3 at room temperature, orthorhombic Fddd at 10 K; supercage Li split over 2-fold-coordinated sites III and III′ | The definitive structure, but paywalled. The ACS supporting-information page returned 403. Its coordinates are in ICSD (subscription), not COD. **Not obtained.** |
| Feuerstein & Lobo, *Chem. Mater.* **10** (1998) 2197 (Li-LSX/LiX, neutron + NMR) | Li-LSX structure | Paywalled. **Not obtained.** |
| COD search (Li–Al–Si–O, 4 and 5 elements) | Only two faujasite entries | See next two rows |
| COD 1530083: Herden, Mortier, Schöllner, Einicke, Uytterhoeven, Gellens, *Zeolites* **2** (1982) 131 | Li-X, Si/Al = 1.23 (Li62.4 H23.6 Al86 Si106 O384), Fd-3m, a = 24.714 Å, Li at SI′ and SII only, Si/Al disordered | Not LSX. Used for starting SI′/SII positions, the starting cell and validation |
| COD 1529943: Forano, Slade, Krogh Andersen, Krogh Andersen, Prince, *J. Solid State Chem.* **82** (1989) 95 | Anhydrous Li-Y (Al57 Si135, Si/Al = 2.37), neutron, Fd-3m, a = 24.450 Å | Not LSX. Used for validation |
| IZA-SC FAU CIF (Baerlocher & McCusker, Database of Zeolite Structures) | Pure-SiO₂ framework, DLS-optimised, Fd-3m origin 2, a = 24.345 Å | Framework source |

All downloaded files are in `sources/`, unmodified. If you have ICSD access,
the Plévert 1997 room-temperature Fd-3 entry should replace this model as the
experimental reference. `03_analyze.py` can compare against it once it is added
to `sources/`.

## 2. How the cell was built (`01_build.py`)

1. **Framework.** The IZA FAU CIF expanded to P1 gives 192 T + 384 O. The cell
   is scaled to a = 24.714 Å (Herden Li-X) as the starting point.
2. **Si/Al = 1 ordering.** The FAU T–T graph is bipartite, so a two-colouring
   gives 96 Si + 96 Al with every O bridging one Si and one Al (Loewenstein's
   rule; asserted in the script). The ordered framework has space group
   **Fd-3 (203)**, the standard description of LSX (checked with spglib).
3. **Li sites, in Fd-3m/Fd-3 origin 2.** Here 16c (0,0,0) is the D6R centre,
   8a (⅛,⅛,⅛) the sodalite centre and 8b (⅜,⅜,⅜) the supercage centre; the
   script checks these from the framework geometry.
   * **SI′**, 32e (x,x,x), x = 0.0392 (Herden): sodalite-cage side of the D6R 6-ring. 32 Li.
   * **SII**, 32e (x,x,x), x = 0.2275 (Herden): supercage side of the sodalite 6-ring. 32 Li.
   * **SIII**, 48f (x,⅛,⅛), on the −4 axis in front of the sodalite 4-ring,
     facing the supercage. The 4-ring centre found from the framework is at
     x = 0.3030; Li starts 2.0 Å beyond it. **32 of the 48 sites are filled.**
   * Occupancies 32/32/32 follow dehydrated Li-LSX:
     - SI′ and SII are full at high Li exchange, and the remaining Li go to
       SIII (LiNaLSX NMR/neutron study,
       [Microporous Mesoporous Mater. 1998, S138718119800211X](https://www.sciencedirect.com/science/article/abs/pii/S138718119800211X)).
     - Plévert 1997 (abstract) puts the supercage Li of room-temperature
       Li-LSX at 2-fold-coordinated sites III/III′.
4. **Choosing 32 of 48 SIII sites.** All 48f sites are equivalent in Fd-3, so
   only the Li–Li arrangement differs between choices. Simulated annealing on
   a Li–Li 1/r sum (40 runs) gave 38 distinct arrangements, nearly degenerate.
   The three lowest all put 4 SIII Li in front of each sodalite cage.

## 3. Relaxation (`02_relax.py`)

* MACE-MP-0 medium + D3(BJ), env `ml` (mace-torch 0.3.17, ASE 3.29), as in
  `hkust1_mace/`. The model ran in **float32**: float64 needed 6.6 GB of the
  8 GB GPU, while float32 needs 3.3 GB and the energy differs by 2×10⁻⁴ eV.
  The run used a 4 GB allocator cap and `odac25/scripts/recon/gpu_watchdog.sh`.
* Stage 1: atoms only, fixed cell, FIRE to fmax 0.05 eV/Å.
* Stage 2 (**main result**): atoms + cubic cell (hydrostatic strain), BFGS to
  fmax 0.01 eV/Å with `exp_cell_factor = 50`, so the residual pressure is
  below about 0.005 GPa.
* Stage 3: free cell from stage 2, done only to measure distortion (below).
* Three SIII arrangements were relaxed. Cubic-cell energies span 10 meV per
  672-atom cell, which is at numerical noise: c1 and c2 are symmetry-equivalent
  and still differ by 3 meV. The lowest, c0, was kept.
* Each relaxation takes about 3 min.

**Why cubic.** The experimental room-temperature cell is cubic because SIII Li
are disordered. With the cell free, one ordered arrangement distorts the cell:
* c0: a, b, c = 24.82, 24.91, 24.57 Å, angles up to 90.7°, −0.35 eV/cell.
* c1/c2: 24.30, 24.90, 25.08 Å (3 % anisotropy).

Since the distortion depends on which 32 of 48 sites are filled, it comes
from ordering SIII in a single cell at 0 K, not from Li-LSX itself. (Plévert
saw a real orthorhombic distortion only at 10 K.) The free-cell structures
are kept in `build/*_freecell.extxyz`.

## 4. Validation (`03_analyze.py`, numbers from `results/geometry.json`)

| Quantity | Model | Experiment / reference |
|---|---|---|
| a (Å) | 24.764 | Li-X (Si/Al 1.23) 24.714 [Herden]; Li-Y (2.37) 24.450 [Forano]. An LSX cell should be the largest of the three. MACE-MP (PBE-trained) typically overestimates cells by about 0.5–1 % |
| Si–O, Al–O mean (Å) | 1.638, 1.758 | typical zeolite values about 1.61, 1.73–1.74; about 1.5 % long (same PBE trend) |
| T–O–T (°) | mean 132.4, range 120.5–150.8 | |
| SI′ Li–O (3-coordinate) (Å) | 1.886 | 1.908 [Herden], 1.911 [Forano] |
| SII Li–O (3-coordinate) (Å) | 1.948 | 2.027 [Herden], 2.038 [Forano] |
| SI′ height above the 6-ring O plane, into the sodalite cage (Å) | +0.30 | −0.06 [Herden], +0.25 [Forano] |
| SII height above the 6-ring O plane, into the supercage (Å) | +0.29 | +0.43 [Herden], +0.16 [Forano] |
| SIII | 1.8 Å outside the 4-ring centre; each Li 0.41 Å (max 0.72) off the −4 axis; 16 two-coordinate + 16 three-coordinate, Li–O 1.96–2.03 Å | Plévert 1997: supercage Li are 2-fold coordinated at III/III′, which matches the off-axis, 2–3-coordinate SIII found here |
| Framework O vs experimental O (same site, after scaling to the model cell) (Å) | 0.07–0.14 [Herden], 0.12–0.19 [Forano] | differences expected from Si/Al = 1 vs 1.23 / 2.37 |
| Shortest Li–Li (Å) | 4.06 | |

The averaged Fd-3 CIF expands back to 576 framework + 32 + 32 + 48 sites,
with composition Li96 Al96 Si96 O384 (pymatgen round trip). Orbit rms spreads
are 0.08–0.18 Å for framework and SI′/SII atoms and 0.50 Å for SIII.

## 5. Comparison with a published Si/Al-ordered structure (`04_compare_reference.py`)

**Reference:** COD 1540962, in `../reference_structure/`.
* Paper: Shepelev, Anderson, Smolin, *Zeolites* **10** (1990) 61,
  doi:10.1016/0144-2449(90)90095-9.
* Material: partially Li-exchanged X, dehydrated at 275 °C:
  Li41.3 Na24.3 [Al92 Si100 O384] (H₂O)8.6. Si/Al = 1.09; Fd-3 origin 2 with
  ordered Si/Al; a = 24.75 Å.
* Method: single-crystal X-ray, according to later citing papers.
* Which state this file is: the same paper's hydrated state is COD 1540961
  (69 H₂O, a = 24.82 Å; kept in `sources/1540961_Shepelev1990_hydrated.cif`).
  That confirms 1540962 is the dehydrated one.

It is the closest open structure to Li-LSX, but it is not Li-LSX: only 66 of
the 92 charge-balancing cations are located, 24 of them are Na, and no Li is
found in the supercage.

**Overlay.** The model (Fd-3 average) was mapped onto the reference with the
Fd-3m operation that best overlays the frameworks. Fd-3 has index 2 in Fd-3m,
so the same structure can be written with Si and Al swapped. Model O2/O3
correspond to reference O3/O2.

| Quantity | Model (Li96, Si/Al 1.00) | Shepelev 1990 (Li41Na24, Si/Al 1.09) |
|---|---|---|
| a (Å) | 24.764 | 24.750 (+0.06 %) |
| Si, Al position offset (Å) | 0.04, 0.04 | – |
| O1, O2, O3, O4 offset (Å) | 0.08, 0.13, 0.11, 0.17 | – |
| Si–O mean (Å) | 1.635 | 1.650 |
| Al–O mean (Å) | 1.753 | 1.690 (Al site has 4 % Si) |
| T–O1–T (°) | 146.5 | 146.6 |
| T–O2–T / T–O3–T (°) | 124.7 / 125.0 | 126.5 / 133.1 |
| T–O4–T (°) | 135.6 | 143.7 |
| Li SI′: x, Li–O, height into sodalite | 0.0481, 1.88 Å, +0.30 Å | 0.053, 1.99 Å, +0.60 Å (25 Li) |
| Li SII: x, Li–O, height into supercage | 0.2246, 1.95 Å, +0.29 Å | 0.210, 2.05 Å, **−0.32 Å** (16 Li; same 6-rings also hold 11 Na at x = 0.25) |
| Supercage (SIII) | 32 Li at 48f, Li–O 1.99 Å (×2) | no Li located; 8.6 Na at 96g (0.42, 0.10, 0.10), 1.3 Å from model SIII; 26 cations unlocated |

**What this shows**
* **Framework: good agreement.** T atoms are within 0.04 Å and O within
  0.17 Å; the cell differs by 0.06 %.
  * The cell match is partly fortuitous: Na enlarges the reference cell, and
    its higher Si/Al and MACE's PBE-type overestimate push in opposite
    directions.
  * The largest framework differences are T–O3–T and T–O4–T, 8° each. Those
    rings border the sites where the reference has Na and water.
  * Al–O differs by 0.06 Å. The reference Al–O (1.69 Å) is short for ordered
    LSX (typically about 1.73 Å), and the model's is about 1.5 % long, so
    the true value probably lies in between.
* **Li SI′: agrees within the experimental spread.** The offset is 0.21 Å.
  The reference heights of 0.60 Å (this X-ray study), −0.06 Å (Herden Li-X,
  neutron) and 0.25 Å (Forano Li-Y, neutron) bracket the model's 0.30 Å.
* **Li SII: disagrees with this reference** (0.63 Å). Here Li sits just
  inside the sodalite cage, but in both neutron Li-faujasites (+0.43, +0.16 Å)
  it sits on the supercage side, like the model. Most likely this is because
  Li and Na share the SII 6-rings in this mixed-cation crystal, and Li is a
  weak X-ray scatterer. It cannot be settled without a pure Li-LSX refinement.
* **SIII: not tested.** No supercage Li is located in this reference.

The Plévert 1997 neutron structure (ICSD) is still the reference needed to
test SII and SIII. `04_compare_reference.py` takes any Fd-3 origin-2 CIF
with labelled Si/Al sites; change `REF` and the cation handling to reuse it.

## 6. Limitations

* Not an experimental structure. Li positions are the MACE minimum starting
  from literature sites; SIII coordinates are model-derived, because no open
  Li-LSX refinement was available.
* SIII ordering is one choice out of many near-degenerate ones (section 3).
  For adsorption simulations, consider averaging over several SIII
  arrangements or letting Li move (e.g. Li as mobile in GCMC/MD).
* 0 K, no thermal motion; MACE-MP-0 bond lengths are about 1.5 % long.
* No partial charges are assigned yet; a force field for N₂ GCMC still has to
  be chosen.

## Reproduce

```bash
conda run -n ml python -I 01_build.py      # ~10 s, build/LiLSX_init_c{0,1,2}.extxyz
# one GPU job at a time; ~3 min per candidate
for c in 0 1 2; do conda run -n ml python -I 02_relax.py $c; done
conda run -n ml python -I 03_analyze.py    # CIFs + results/
conda run -n ml python -I 04_compare_reference.py   # vs COD 1540962 (Shepelev 1990)
```
