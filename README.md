# SiB4: Carbon Isotopes and Fire Disturbance Development

**Research software | Active development**

Repository: https://github.com/aleya-kaushik/c13c14_fireba_r

## Overview

This repository contains a research-development branch of the Simple Biosphere Model version 4 (SiB4), a process-based terrestrial biosphere model that simulates exchanges of carbon, water, and energy between land ecosystems and the atmosphere. SiB4 represents vegetation growth and seasonal phenology, photosynthesis, respiration, carbon allocation, and carbon storage in living vegetation and dead organic matter.

The development work represented here focuses on carbon-isotope tracking (carbon-13 and carbon-14) and fire-disturbance processes. These capabilities support investigations of how terrestrial carbon is taken up, redistributed, released, and retained over timescales ranging from seasonal vegetation growth to longer-term ecosystem carbon turnover.

## Scientific scope

Current research and development in this code lineage includes:

- **Carbon isotopes:** Representation and diagnostic tracking of carbon-13 and carbon-14 alongside the model's carbon pools and fluxes, for studies of carbon cycling and turnover.
- **Fire disturbance:** Fire-related carbon losses and subsequent redistribution of carbon among vegetation, litter, and soil pools.
- **Alternative fire representations:** An established emissions-driven treatment and a burned-area-driven treatment under development and validation. The latter is intended to distinguish direct combustion from vegetation mortality and transfers into dead carbon pools.
- **Ecosystem recovery:** Evaluation of vegetation regrowth, phenology, photosynthetic activity, and the consequences of disturbance for ecosystem carbon stocks.

The presence of code for a process or diagnostic does not necessarily mean that it is enabled in every run configuration or has completed independent validation.

## Development and validation status

**This is an evolving scientific codebase, not a frozen or fully validated software release.** Model behavior depends on the selected source revision, namelist settings, input data, and simulation protocol.

In particular:

- The burned-area-driven fire treatment is developmental and should not be assumed to be interchangeable with the established emissions-driven configuration.
- Carbon and isotope diagnostics, including derived isotope ratios, should be independently checked against the underlying pool and flux quantities before quantitative interpretation.
- Reproducing a scientific result requires recording the exact commit, compiler and library environment, forcing datasets, parameter files, initialization and spin-up procedure, namelist settings, and output selection.

Development-stage features should not be interpreted as validated scientific findings without accompanying tests and documentation.

## Obtaining the code

Clone the repository using Git:

`git clone https://github.com/aleya-kaushik/c13c14_fireba_r.git`

`cd c13c14_fireba_r`

For reproducible use, record the exact commit with `git rev-parse HEAD`.

**Build and run instructions:** SiB4 is a compiled Fortran model with external input data and configuration requirements. The precise build commands, compiler dependencies, input-file layout, and example run settings for this published repository still need to be documented and verified against its source tree. Please do not assume that cloning the repository alone provides a complete, runnable experiment.

## Using and interpreting this repository

This repository can document ongoing model development and serve as a reference for specific code implementations. Scientific results should identify which features were active, distinguish established processes from experimental alternatives, and report relevant validation checks. Results obtained with different fire schemes or initialization protocols should not be treated as directly comparable without documenting those differences.

## Citation and attribution

Until a formal archived software release is available, cite the repository using its URL and a specific commit hash or tagged release, together with the appropriate original SiB4 publications.

An example citation format is:

SiB4 carbon-isotope and fire-disturbance development code. GitHub repository, https://github.com/aleya-kaushik/c13c14_fireba_r (accessed YYYY-MM-DD; commit: [hash]).

This repository builds on SiB4, including the foundational model descriptions:

- Haynes, K. D., et al. (2019). Representing grasslands using dynamic prognostic phenology based on biological growth stages: 1. Implementation in the Simple Biosphere Model (SiB4). *Journal of Advances in Modeling Earth Systems*. https://doi.org/10.1029/2018MS001540
- Haynes, K. D., et al. (2019). Representing grasslands using dynamic prognostic phenology based on biological growth stages: Part 2. Carbon cycling. *Journal of Advances in Modeling Earth Systems*. https://doi.org/10.1029/2018MS001541

Please preserve attribution to the original SiB4 developers and contributors when discussing or reusing this research branch.

## Licensing and reuse

Public visibility on GitHub does not, by itself, establish an open-source license. Consult any LICENSE file and relevant upstream SiB4 terms before reusing, redistributing, or incorporating the code into another software package. Where no applicable license is supplied, contact the repository maintainer about permissions.

## Contact

For questions or discussion, use the repository's GitHub Issues feature if enabled, or contact the repository maintainer through their GitHub profile.
