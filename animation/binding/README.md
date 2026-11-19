# Phase 2: N₂ end-on at an exposed SIII Li⁺ (binding geometry for the animation)

**What this is.** It gives a plausible bound geometry for the scripted animation,
nothing more. The N₂ was placed end-on at one SIII Li⁺, then only the two N atoms
were relaxed with MACE-MP-0 medium + D3(BJ) (float32, env `ml`). The framework
and all other Li stayed fixed.

## Choice of Li

The supercage is centred at (⅜,⅜,⅜), the same cage rendered in Phase 1. Its
four SIII Li are atoms 662, 654, 652 and 651 (0-based, in
`structure/LiLSX_P1_MACE.cif`). **Li 662** was chosen because:
* it is 2-fold O-coordinated (Li–O 1.92 and 1.96 Å; the next O is at 2.84 Å),
  so it is the most exposed;
* it sits on the back half of the cage, so the Phase-4 cutaway keeps it.

## Result (`binding_summary.json`)

| Quantity | MACE-MP-0 + D3, framework fixed |
|---|---|
| Li–N | **2.063 Å** |
| Li–N–N | **171.2°**, end-on and nearly linear; the N₂ leans toward the wall |
| Li to N₂ centre of mass | 2.616 Å |
| N–N | 1.1158 Å (MACE gas-phase N₂: 1.1124 Å) |
| Closest other contacts | O at 3.25–3.38 Å; the outer N is 3.51 Å from a neighbouring SII Li (617) |
| Interaction energy vs. gas N₂ | −33.7 kJ/mol (MACE; see the comparison below) |

**Is the result independent of the starting geometry?**
* Linear starts at Li–N 2.2, 2.6 and 3.0 Å all relax to 2.02 Å and 179°.
  Symmetry cancels the sideways force on a perfectly linear N₂, so these runs
  stay linear.
* A scan of 16 tilted starts (`02_tilt_scan.py`: 15° and 30° tilts × 8
  directions) finds two bent basins, at 171° (the minimum) and 169° (+3.3 kJ/mol).
* The exactly linear geometry is **8.8 kJ/mol above** the minimum, so it is a
  ridge, not a minimum. The lean puts the outer N near the neighbouring SII Li
  and the wall O, which MACE + D3 rewards.

## Comparison with published values

* **Gas-phase Li⁺–N₂ ab initio** (bare cation, linear minimum):
  * 2018 MRCI/aug-cc-pVQZ potential surface (*J. Chem. Sci.* 130, 155,
    doi:10.1007/s12039-018-1561-x): global minimum collinear with
    R = 4.96 a₀ = **2.62 Å**. A 1975 SCF+CI surface (*Chem. Phys.*,
    ScienceDirect 030101047585021X) has R ≈ 5.0 a₀ (2.65 Å) from Li⁺ to the N₂
    centre of mass.
  * If R is the Li⁺–centre-of-mass distance (the usual Jacobi coordinate), the
    Li–N contact is R − r/2 ≈ 2.62 − 0.55 = **≈ 2.07 Å**. The model's Li to
    centre of mass (2.616 Å) and Li–N (2.063 Å) match this almost exactly.
  * Caveat: both values come from abstracts and search snippets. The full
    texts were behind logins, so the definition of R is not confirmed.
  * MP2/CCSD on Li⁺(N₂)ₙ (*Struct. Chem.* 2018, doi:10.1007/s11224-018-1147-8):
    the n = 1 complex is linear, with a binding energy of 11.1 kcal/mol
    (46 kJ/mol).
* **Li-LSX experiment:** the low-coverage N₂ heat of adsorption is ≈ 22.6 kJ/mol
  (LiNa-LSX heats study, *Chem. Eng. Sci.* 2017, ScienceDirect
  S0009250917307765). MACE + D3 gives 33.7 kJ/mol, i.e. it **overbinds by
  about 50 %**.
* **No published Li-LSX value found.** No published framework (zeolite) Li–N
  distance turned up. Kowsari 2018 (*Microporous Mesoporous Mater.* 264), which
  is the likely source, studies N₂ at Li-III in Li-LSX but its numbers weren't
  accessible.

**Reading.** A framework Li⁺ is partly screened by its O neighbours, so it
should hold N₂ a little more loosely than a bare Li⁺. The expected Li–N is then
the gas-phase ≈ 2.07 Å or somewhat longer, roughly 2.1–2.3 Å. MACE's 2.06 Å sits
at the short end, consistent with its overbinding.

The 2.2–2.6 Å guess spans both conventions: 2.6 Å is the Li-to-N₂-centre
distance, while the Li–N contact is about 0.55 Å shorter. For the animation:
* 2.06 Å and 171° are a plausible docked geometry;
* the 2.8 Å glow threshold safely brackets it.

## Files

* `01_bind_n2.py`: Li choice, clearance axis, 4 starts, rough energies → `binding.json`, `starts/`
* `02_tilt_scan.py`: 16 tilted starts → `tilt_scan.json`; the minimum goes to `complex_relaxed.extxyz`
* `03_summary_view.py` (env `viz`): `binding_summary.json` and `bound_view.png` (Phase-1 cutaway camera, active Li highlighted)
