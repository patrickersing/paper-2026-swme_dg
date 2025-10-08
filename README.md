# paper-2025-swme-dg_dev

Reproducibility repository for the upcoming paper "Entropy stable discontinuous Galerkin schemes
for the shallow water moment equations"

## Installation
1. Install Julia v1.10
2. Download the repository
```
git clone https://github.com/patrickersing/paper-2025-swme-dg_dev.git
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

Then execute the following commands to run e.g. "elixir_shallowwater_linearized_moments_shock_capturing.jl":
```julia
include("examples/elixir_shallowwater_linearized_moments_shock_capturing.jl")
```
To visualize the results with Plots.jl run
```julia
using Plots
plot(sol)
```