
using OrdinaryDiffEqSSPRK, OrdinaryDiffEqLowStorageRK
using Trixi
using TrixiShallowWater

include("../src/main.jl")

# Semidiscretization of the shallow water moment equations
equations = ShallowWaterLinearizedMomentEquations1D(gravity = 9.812, H0 = 1.75,
                                                    n_moments = 2)

function initial_condition_stone_throw(x, t,
                                       equations::Union{ShallowWaterMomentEquations1D,
                                                        ShallowWaterLinearizedMomentEquations1D})
    # Initial lake-at-rest configuration
    H = 1.75
    v = 0.0
    a = zeros(equations.n_moments)

    # Set discontinuous velocity / moments
    eps = 1e-3
    if x[1] >= -1.0 && x[1] <= 0.0
        v = eps
        a = eps * ones(equations.n_moments)
    elseif x[1] > 0.0 && x[1] <= 1.0
        v = -eps
        a = -eps * ones(equations.n_moments)
    end

    # Set smooth bottom topography
    b = 1 / exp(0.5 * x[1]^2)

    h = H - b

    return SVector(h, h * v, (h * a)..., b)
end

initial_condition = initial_condition_stone_throw

###############################################################################
# Get the DG approximation space

volume_flux = (flux_ec, flux_nonconservative_ec)
surface_flux = (FluxPlusDissipation(flux_ec,
                                    DissipationLaxFriedrichsEntropyVariables(Trixi.max_abs_speed)),
                flux_nonconservative_ec)

indicator_var(u, equations) = u[2]^3
basis = LobattoLegendreBasis(1)
indicator_sc = IndicatorHennemannGassner(equations, basis,
                                         alpha_max = 0.5,
                                         alpha_min = 0.001,
                                         alpha_smooth = true,
                                         variable = indicator_var)

volume_integral = VolumeIntegralShockCapturingHG(indicator_sc;
                                                 volume_flux_dg = volume_flux,
                                                 volume_flux_fv = surface_flux,)

solver = DGSEM(basis, surface_flux, volume_integral)

###############################################################################
# Create the TreeMesh for the domain [-1, 1]

coordinates_min = -4.0
coordinates_max = 4.0

mesh = TreeMesh(coordinates_min,
                coordinates_max,
                initial_refinement_level = 6, # 2^refinement_level
                n_cells_max = 10_000,
                periodicity = true)

# create the semi discretization object
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver)

###############################################################################
# ODE solver
tspan = (0.0, 8000.0)
ode = semidiscretize(semi, tspan)

###############################################################################
# Callbacks
summary_callback = SummaryCallback()

analysis_interval = 10000
analysis_callback = AnalysisCallback(semi, interval = analysis_interval,
                                     save_analysis = true,
                                     extra_analysis_integrals = (entropy,
                                                                 lake_at_rest_error))
alive_callback = AliveCallback(analysis_interval = analysis_interval)
save_solution = SaveSolutionCallback(dt = 100.0, save_initial_solution = true,
                                     save_final_solution = true)
stepsize_callback = StepsizeCallback(cfl = 0.9)

callbacks = CallbackSet(summary_callback, analysis_callback, alive_callback, save_solution,
                        stepsize_callback)
###############################################################################
# run the simulation

sol = solve(ode,
            CarpenterKennedy2N54(williamson_condition = false);
            dt = 1.0,              # solve needs some value here but it will be overwritten by the stepsize_callback
            ode_default_options()...,
            callback = callbacks,
            maxiters = 10_000_000,);
