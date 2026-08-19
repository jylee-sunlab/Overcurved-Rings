# Planar instability and reconfiguration of overcurved rings

## Files

- `run_overcurved_ring.m` — main script. Edit the input block and run the complete analysis.
- `overcurved_energy_scan.m` — performs the nonlinear equilibrium scan over the prescribed overcurving ratio.
- `save_scan_results.m` — writes the scan results and optional field data to CSV or XLSX files.
- `plot_scan_results.m` — plots the realized overcurving ratio and energy partition.
- `plot_shapes.m` — reconstructs and saves three-dimensional ring shapes and the progression panel.
- `reference/` — reference CSV outputs for supplied benchmark cases.
- `CITATION.cff` — citation metadata.
- `LICENSE` — MIT license.

## Usage

The simplest workflow is to edit the input section of `run_overcurved_ring.m` and execute the script.

```matlab
run_overcurved_ring
```

## Model

The prescribed overcurving ratio is

```text
O_p = m * theta_0 / pi
```

where `theta_0 = kappa_p * L0` is the intrinsic angular span of one undeformed lobe.

For the present implementation, `m = 2`, so

```text
1 <= O_p <= 3
```

The equilibrium ring is obtained by minimizing the elastic energy over the global variables `omega` and `ell` and the arc-length fields `theta_m(s0)` and `q_map(s0)`.

The coded total-ring energy is equivalent to

```text
U = m * integral_0^L0 [
      E*I1*(kappa_m1 - kappa_p)^2
    + E*I2*kappa_m2^2
    + G*J*tau_m^2
    + E*A*eps_0^2
] ds0
```

with

```text
G = E / (2*(1 + nu)).
```

Here

- `kappa_m1` and `kappa_m2` are the material bending curvatures.
- `tau_m` is the mechanical torsion.
- `eps_0` is the axial strain.
- `theta_m(s0)` is the material-frame rotation.
- `q_map(s0)` generates the monotone material mapping `t(s0)`.
- `ell` is the geometric scale factor.
- `omega` is the lobe angular span.

The principal output is the realized overcurving ratio `O_geom`, stored as `Oeq` in the result table.

## Input parameters

All dimensional quantities use SI units.

| Field | Meaning |
| --- | --- |
| `m` | lobe-pair count. The present implementation requires `m = 2` |
| `E` | Young's modulus |
| `nu` | Poisson's ratio |
| `I1` | second moment of area of the bending channel carrying the preset curvature |
| `I2` | second moment of area of the second bending channel |
| `J` | Saint-Venant torsional constant |
| `A` | cross-sectional area |
| `kp` | preset curvature `kappa_p` |
| `Ns0` | spatial grid points per representative lobe. Must be odd |
| `nDivTheta0` | number of scan intervals over `O_p` |

The default numerical settings in `overcurved_energy_scan.m` may be overridden through additional fields of `p`.

## Output

`overcurved_energy_scan` returns a structure `res`.

### `res.table`

One row is stored for each prescribed overcurving ratio.

| Column | Meaning |
| --- | --- |
| `Op` | prescribed overcurving ratio |
| `theta0` | undeformed lobe angular span |
| `L0` | undeformed lobe length |
| `omega` | equilibrium lobe angular span |
| `phi` | geometric angle of the lobe construction |
| `ell` | equilibrium scale factor |
| `Oeq` | realized overcurving ratio `O_geom` |
| `Ub` | bending energy |
| `Ut` | torsional energy |
| `Ua` | axial energy |
| `Utotal` | total ring energy |
| `exitflag` | `fmincon` convergence flag |

### `res.fields`

The converged arc-length fields are stored for every scan point.

- `s0`
- `theta_m`
- `t`
- `dt_ds0`
- `lambda`
- `q`

The result structure also contains solver diagnostics in `res.diag`, validated parameters and derived rigidities in `res.p`, and run metadata in `res.meta`.

## Ring shapes

`plot_shapes` reconstructs the closed three-dimensional centerline associated with a prescribed realized ratio `O_geom`.

Each ring is assembled from symmetry-related constant-curvature lobe segments.
The shape construction is restricted to `m = 2`.

The endpoints

```text
O_geom = 1
O_geom = 3
```

are singular in the closed-form shape parameterization.
The default shape range therefore stops short of both endpoints.

The reconstruction also reports a closure residual.
For a valid shape this residual should remain near numerical round-off.

## Requirements

- Developed and tested with MATLAB R2025b.
- MATLAB Optimization Toolbox is used.
- `fmincon` is used for the equilibrium minimization.
- `fsolve` is used by `plot_shapes`.


## Numerical notes

The solver periodically probes the opposite energy basin and retains the lower-energy solution.

For quantitative calculations, `Ns0 = 201` or finer is recommended.

For non-circular cross-sections, use the Saint-Venant torsional constant `J` rather than the polar second moment.


## Reference data

The `reference/` folder contains numerical outputs for supplied benchmark cases.

- `rect_4x3_kp16_Ns101_summary.csv`
- `rect_4x3_kp16_Ns101_parameters.csv`
- `ellipse_12x3_kp16_Ns201_summary.csv`

These files can be used to check the output of a local installation against previously generated results.


## Citation

If you use this code, please cite:

> Jihoon Yook and Jae Young Lee, “Planar instability and reconfiguration of overcurved rings,”
> *Submitted*

## License

MIT License — see [LICENSE](LICENSE).
 