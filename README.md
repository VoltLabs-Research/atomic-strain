# AtomicStrain

`AtomicStrain` computes atomic strain metrics from a current frame and an optional reference frame.

## CLI

Usage:

```bash
atomic-strain <lammps_file> [output_base] [options]
```

### Arguments

| Argument | Required | Description | Default |
| --- | --- | --- | --- |
| `<lammps_file>` | Yes | Input LAMMPS dump file. | |
| `[output_base]` | No | Base path for output files. | derived from input |
| `--cutoff <float>` | No | Cutoff radius for neighbor search. | `3.0` |
| `--reference <file>` | No | Reference LAMMPS dump file. If omitted, the current frame is used. | current frame |
| `--eliminateCellDeformation` | No | Eliminate cell deformation before computing strain. | `false` |
| `--assumeUnwrapped` | No | Assume coordinates are already unwrapped. | `false` |
| `--calcDeformationGradient` | No | Compute deformation gradient `F`. | `true` |
| `--calcStrainTensors` | No | Compute strain tensors. | `true` |
| `--calcD2min` | No | Compute `D²min`. | `true` |
| `--threads <int>` | No | Maximum worker threads. | auto |
| `--help` | No | Print CLI help. | |

## Build With CoreToolkit

```bash
cd /path/to/voltlabs-ecosystem/tools/CoreToolkit
conan create . -nr

cd /path/to/voltlabs-ecosystem/plugins/AtomicStrain
conan create . -nr
```
