# AtomicStrain

`AtomicStrain` computes atomic strain metrics from a current frame and an optional reference frame.

## One-Command Install

```bash
curl -sSL https://raw.githubusercontent.com/VoltLabs-Research/CoreToolkit/main/scripts/install-plugin.sh | bash -s -- AtomicStrain
```

## Build from source

Requires [Conan 2.x](https://docs.conan.io/2/installation.html), CMake 3.20+, and a C++23 compiler (GCC 14+ or Clang 17+).

### Prerequisites

The following Conan packages must be available in your local cache:

- `coretoolkit/1.0.0` (from the `CoreToolkit` repository)

For each dependency, clone its repository and create the package:

```bash
conan create <path-to-dependency-repo> --build=missing -o "hwloc/*:shared=True"
```

### Build

From the root of this repository:

```bash
conan install . -of build --build=missing -o "hwloc/*:shared=True"
cmake --preset conan-release
cmake --build build/build/Release -j
```

### Run

```bash
./build/build/Release/atomic-strain --help
```

### Package as Conan recipe

To make this plugin available as a Conan package for other projects:

```bash
conan create . --build=missing -o "hwloc/*:shared=True"
```

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
| `--eliminate_cell_deformation` | No | Eliminate cell deformation before computing strain. | `false` |
| `--assume_unwrapped` | No | Assume coordinates are already unwrapped. | `false` |
| `--calc_deformation_gradient` | No | Compute deformation gradient `F`. | `true` |
| `--calc_strain_tensors` | No | Compute strain tensors. | `true` |
| `--calc_d2min` | No | Compute `D²min`. | `true` |
| `--threads <int>` | No | Maximum worker threads. | auto |
| `--help` | No | Print CLI help. | |
