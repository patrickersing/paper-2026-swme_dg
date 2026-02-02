# Entropy analysis and entropy stable DG methods for the Shallow Water Moment Equations

[![License: MIT](https://img.shields.io/badge/License-MIT-success.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/1043271198.svg)](https://doi.org/10.5281/zenodo.18460096)

This reproducibility repository contains the information and code to reproduce the results of the article "Entropy analysis and entropy stable DG methods for the Shallow Water Moment Equations"

If you find these results useful, please cite the article mentioned above. If you use the implementations provided here, please also cite this repository as

```bibtex
@software{ersing2026swmeRepro,
  title={Reproducibility repository for "Entropy analysis and entropy stable {DG} methods for the {S}hallow {W}ater {M}oment {E}quations"},
  author={Julio Careaga, Patrick Ersing},
  year={2026},
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.18460096},
  url          = {https://doi.org/10.5281/zenodo.18460096},
}
```

## Installation
1. Install Julia v1.10
2. Download the repository
```
git clone git@github.com:patrickersing/paper-2025-swme-dg_dev.git
```
1. Set the working directory to the main folder and instantiate the Julia environment using
```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

## Usage
To reproduce the results in the paper start a new `Julia` session with the `--project=.` flag. It is possible (but optional) to use multiple threads by appending the `--threads=N` flag, to run on `N` threads. 
```bash
julia --project=. --threads=N
```

Then execute the following commands in the Julia-REPL to create the respective results:

- **Figure 1**
The following command will create Figure 1 and save it as a `.pdf`
```julia 
include("/figures/plot_figure_1.jl")
```

- **Figure 2**
The following command will create Figure 2 and save it as a `.pdf`
```julia 
include("/figures/plot_figure_2.jl")
```

- **Figure 3**
The following command will create Figure 3 and save it as a `.pdf`
```julia 
include("/figures/plot_figure_3.jl")
```

- **Figure 4 & 5**
The following command will create Figures 4 & 5 and and save them as a `.pdf`
```julia 
include("/figures/plot_figure_4.jl")
```

- **Table 1 & 2**
The following command will print the convergence results for Table 1 to the Julia-REPL
```julia
using Trixi
convergence_test("/examples/elixir_shallowwater_moments_convergence.jl", 5)
convergence_test("/examples/elixir_shallowwater_linearized_moments_convergence.jl", 5)
```

## Authors
- [Patrick Ersing](https://liu.se/en/employee/pater53) (Linköping University, Sweden)
- [Julio Cesar Careaga Solis](https://www.rug.nl/staff/j.c.careaga.solis/?lang=en) (University of Groningen, Netherlands)

## Disclaimer

Everything is provided as is and without warranty. Use at your own risk!
