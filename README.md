# Atomic Strain

Computes per-atom strain (deformation gradient, strain tensors, D²min) relative to a reference configuration.

## Install

```bash
vpm install @voltlabs/atomic-strain
```

## CLI

```bash
atomic-strain <input_dump> [output_base] [options]
```

| Argument | Required | Default | Description |
|---|---|---|---|
| `<input_dump>` | yes | — | Input LAMMPS dump. |
| `[output_base]` | no | derived from input | Base path for output files. |
| `--reference <file>` | no | current frame | Reference LAMMPS dump frame; if omitted, the current frame is used. |
| `--cutoff <float>` | no | `3.0` | Cutoff radius for neighbor search. |
| `--eliminate_cell_deformation` | no | `false` | Eliminate cell deformation before computing strain. |
| `--assume_unwrapped` | no | `false` | Assume coordinates are already unwrapped. |
| `--calc_deformation_gradient` | no | `true` | Compute the deformation gradient `F`. |
| `--calc_strain_tensors` | no | `true` | Compute strain tensors. |
| `--calc_d2min` | no | `true` | Compute `D²min`. |
| `--calc_polar_decomp` | no | `true` | Compute the polar decomposition (rotation `R` and stretch `U` tensors). |
| `--threads <int>` | no | auto | Maximum worker threads. |

## Exports

| Output file | Exposure | Exporter → artifact |
|---|---|---|
| `{output_base}_atomic_strain.parquet` | Atomic Strain | — (listing-only) |
| `{output_base}_atoms.parquet` | Atomic Strain Model | AtomisticExporter → glb |

---

Full input contract and examples: https://docs.voltcloud.dev/docs/plugins/atomic-strain
