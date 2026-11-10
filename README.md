# CCpol23+(2025) MBMD

Reference implementation of the CCpol23+(2025) water potential energy surface (PES), with a Fortran/MPI molecular dynamics implementation and parameter sets for two three-body representations.

The model updates the three-body interaction while retaining the remaining potential components from the framework introduced in:

> Góra et al., “Predictions for water clusters from a first-principles two- and three-body force field,” The Journal of Chemical Physics 140, 194101 (2014).
> DOI: 10.1063/1.4875097 (https://doi.org/10.1063/1.4875097)

## Related Publication

The optimization of the three-body representation and the analysis
of long-range three-body interactions are described in:

> **Ommair Ishaque, John W. Melkumov, and Krzysztof Szalewicz**,  
> “How Long-Range Are Three-Body ‘Exchange’ Interactions in Liquid Water?”  
> *Journal of Chemical Theory and Computation* **22** (6), 2757–2781 (2026).  
> DOI: [10.1021/acs.jctc.5c02107](https://doi.org/10.1021/acs.jctc.5c02107)

## Overview

CCpol23+(2025) provides an improved representation of short- and long-range three-body interactions in water. Separate short-range (SR) and long-range (LR) fits are connected through a smooth switching function.

Two fitted three-body variants are included:

| Variant | Description | Input directory |
|---|---|---|
| FS3 | Improved $F_{S³}$ representation | input_fit_ccpol23_2025_fs3 (input_fit_ccpol23_2025_fs3) |
| FS2+FS3 | Optimized combination of $F_{S²}$ and $F_{S³}$ | input_fit_ccpol23_2025_fs2fs3 (input_fit_ccpol23_2025_fs2fs3) |

Here, $F_{S²}$ and $F_{S³}$ denote contributions to the three-body interaction. They should not be confused with the two-body and three-body terms of the many-body expansion.

The input directories contain simulation inputs and fitted parameter sets. The molecular dynamics implementation uses rigid monomers, with molecular translation and rotation.

## Short- and Long-Range Representations

The switching coordinate is the average separation between the centers of mass of the three monomers:

$R_{avg}$ = $(R_{AB}+R_{BC}+R_{CA})$/3

The SR and LR fitting domains overlap:

- SR fitting domain: $R_{avg}$ < 5.5 Å.
- LR fitting domain: $R_{avg}$ > 4.5 Å.

During evaluation, the FS3 contribution is selected as follows:

| Average monomer separation | FS3 evaluation |
|---|---|
| $R_{avg}$ ≤ 4.5 Å | SR representation |
| 4.5 < $R_{avg}$ < 5.5 Å | Smooth combination of SR and LR representations |
| $R_{avg}$ ≥ 5.5 Å | LR representation |

## Repository Structure

CCpol23_plus_2025_MBMD/
├── README.md
├── source/
│   ├── Makefile
│   ├── md.f90
│   ├── initialization.f90
│   ├── common_arrays.f90
│   ├── threebody_potential.f90
│   ├── nb_induction_model.f90
│   ├── ewald.f90
│   ├── evolve.f90
│   ├── statistics.f90
│   ├── utility.f90
│   ├── molecular_sites.f90
│   ├── vdw_*.f90
│   └── tests/
│       └── check_fs2_periodic.py
├── input_fit_ccpol23_2025_fs3/
│   ├── CONTROL
│   ├── FIELD
│   ├── Random.num
│   ├── data3b_full
│   ├── data3b_full_SR_5p5
│   ├── data3b_full_LR_4p5
│   └── run.sh
└── input_fit_ccpol23_2025_fs2fs3/
    ├── CONTROL
    ├── FIELD
    ├── Random.num
    ├── data3b_full
    ├── data3b_full_SR_5p5
    ├── data3b_full_LR_4p5
    └── run.sh

## Requirements

The molecular dynamics executable requires:

- A Fortran compiler with compatible MPI Fortran bindings.
- An MPI implementation providing mpifort and mpirun.
- make.
- A shell environment for the example commands.

The supplied job scripts use the Slurm workload manager.

The optional FS2 validation test additionally requires:

- Python 3.
- GNU Fortran (gfortran).

Example instructions for installing Open MPI 4.1.8 are provided in source/readme (source/readme). Use compatible MPI installations when compiling and running the executable.

## Compilation

Clone the repository and build the executable:

git clone https://github.com/ommair/CCpol23_plus_2025_MBMD.git
cd CCpol23_plus_2025_MBMD

make -C source FC=mpifort

The resulting executable is:

source/execute

To rebuild all source files:

make -B -C source FC=mpifort

## Parameter Files

The current parameter reader uses the following filenames in the simulation working directory:

| File | Role |
|---|---|
| CONTROL | Simulation settings, interaction switches, cutoffs, and restart options |
| FIELD | Molecular geometry, interaction-site definitions, and model parameters |
| Random.num | Random-number data supplied with the examples |
| data3b_full | Active SR FS3 parameters; also initializes legacy FS3 arrays |
| data3b_full_LR_4p5 | Active LR FS3 parameters |
| REVCON | Restart state, required when restarting is enabled |

The current reader does not automatically load data3b_full_SR_5p5. To use that named SR parameter set, copy it to the active filename data3b_full in the run directory.

Some FS2 constants are defined directly in source/threebody_potential.f90. Changing input files alone therefore does not replace every FS2 parameter.

Keep the SR and LR parameter files from the same model variant together.

## Model Selection

The three-body control line contains three logical values in this order:

nonadditive_3B   S3_term   S2_term

The supplied variants use:

| Variant | Control values |
|---|---|
| FS3 | T T F |
| FS2+FS3 | T T T |

CONTROL is a positional input file. Preserve the order of its records when editing simulation settings.

## Quick Start: FS3 Variant

Run simulations in a separate working directory so that the supplied examples remain available for reference.

From the repository root:

mkdir -p runs/fs3_example

cp input_fit_ccpol23_2025_fs3/CONTROL \
   input_fit_ccpol23_2025_fs3/FIELD \
   input_fit_ccpol23_2025_fs3/Random.num \
   input_fit_ccpol23_2025_fs3/data3b_full_LR_4p5 \
   runs/fs3_example/

cp input_fit_ccpol23_2025_fs3/data3b_full_SR_5p5 \
   runs/fs3_example/data3b_full

This explicitly stages the named SR fit under the filename expected by the current reader.

Before launching, edit runs/fs3_example/CONTROL:

- Set the restart flag to F for a fresh simulation.
- Retain T T F for the three-body switches.
- Choose the desired temperature, pressure, timestep, and number of steps.
- Use a small step count for an initial execution check.
- Select an appropriate equilibration period for a production simulation.

The supplied FS3 example currently enables restarting. A fresh run requires changing that flag because a REVCON restart file is not included.

Launch a single-rank MPI calculation:

cd runs/fs3_example
mpirun -np 1 ../../source/execute

Increase the MPI rank count as appropriate for the available resources.

## Running the FS2+FS3 Variant

From the repository root:

mkdir -p runs/fs2fs3_example

cp input_fit_ccpol23_2025_fs2fs3/CONTROL \
   input_fit_ccpol23_2025_fs2fs3/FIELD \
   input_fit_ccpol23_2025_fs2fs3/Random.num \
   input_fit_ccpol23_2025_fs2fs3/data3b_full_LR_4p5 \
   runs/fs2fs3_example/

cp input_fit_ccpol23_2025_fs2fs3/data3b_full_SR_5p5 \
   runs/fs2fs3_example/data3b_full

Check the simulation settings and retain T T T for the three-body switches. For a fresh simulation, the restart flag must be F.

Then run:

cd runs/fs2fs3_example
mpirun -np 1 ../../source/execute

## Input Units

The main simulation inputs use the following conventions:

| Quantity | Unit |
|---|---|
| Temperature | K |
| Timestep | ps |
| Two-body and three-body cutoff distances | Å |
| Pressure entered in CONTROL | katm |

For example:

0.001

in the pressure record corresponds to 1 atm.

The pressure comments in the supplied CONTROL examples currently say atm; the implementation reads this value in katm and converts it internally.

Parameter-file coefficients follow the conventions of their corresponding Fortran routines. Their units should be checked against those definitions before modifying or converting fitted parameters.

## Slurm Submission

The supplied run.sh files are cluster-specific examples. Adapt the partition, MPI installation, resource requests, and executable path to the target system.

The current scripts contain hard-coded executable paths. To make the executable configurable, replace the active assignment with:

EXECUTABLE="${EXECUTABLE:-./execute}"

With that assignment in place, an executable path can be supplied at submission:

sbatch --export=ALL,EXECUTABLE=/absolute/path/to/source/execute run.sh

Submit from the simulation working directory containing the required input and parameter files.

## Output Files

The molecular dynamics program writes:

| File | Description |
|---|---|
| OUTPUT | Simulation settings and reported thermodynamic quantities |
| HISTORY.xyz | Saved molecular configurations |
| DIPOLE | Dipole output |
| REVCON | State used to continue a simulation |

Additional diagnostic files may be produced when force or pressure checks are enabled.

Refer to the output headers and source definitions for the units and normalization of individual columns. Distinguish simulation-cell quantities from quantities normalized per molecule when performing subsequent analysis.

Choose production lengths, equilibration periods, and statistical block sizes according to the observable being calculated. The small step counts in the supplied input examples are demonstration settings.

## Restarting a Simulation

To continue an existing simulation:

1. Retain the generated REVCON file in the working directory.
2. Keep the corresponding model parameters and molecular definitions.
3. Enable the restart flag in CONTROL.
4. Check the step limit and other continuation settings before launching.

The implementation appends to several output files. Use separate working directories for independent simulations.

## Citation

> Ishaque, O.; Melkumov, J. W.; Szalewicz, K.  
> How Long-Range Are Three-Body “Exchange” Interactions in Liquid Water?  
> *J. Chem. Theory Comput.* **2026**, *22* (6), 2757–2781.  
> DOI: [10.1021/acs.jctc.5c02107](https://doi.org/10.1021/acs.jctc.5c02107)

