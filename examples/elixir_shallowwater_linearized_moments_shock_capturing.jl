
using OrdinaryDiffEqSSPRK, OrdinaryDiffEqLowStorageRK
using Trixi
using TrixiShallowWater

include("../src/main.jl")

###############################################################################
# Semidiscretization of the shallow water linearized moment equations

equations =
    ShallowWaterLinearizedMomentEquations1D(gravity = 9.812, H0 = 1.75, n_moments = 4)

# Initial condition with a truly discontinuous velocity and smooth bottom topography.
function initial_condition_stone_throw_discontinuous_bottom(
    x,
    t,
    equations::ShallowWaterLinearizedMomentEquations1D,
)

    # Calculate primitive variables

    # flat lake
    H = equations.H0

    # Discontinuous velocity
    v = 0.0
    if x[1] >= -0.75 && x[1] <= 0.0
        v = -0.3
    elseif x[1] >= 0.0 && x[1] <= 0.75
        v = 0.3
    end

    a = ones(equations.n_moments)

    b = cos(0.5 * π * x[1]) / 10.0

    return prim2cons(SVector(H, v, a..., b), equations)
end

initial_condition = initial_condition_stone_throw_discontinuous_bottom

###############################################################################
# Get the DG approximation space

volume_flux = (flux_ec, flux_nonconservative_ec)
surface_flux = (flux_ec, flux_nonconservative_ec)
#surface_flux = (FluxPlusDissipation(flux_ec, DissipationLocalLaxFriedrichs(Trixi.max_abs_speed)), flux_nonconservative_ec)

solver = DGSEM(
    polydeg = 3,
    surface_flux = surface_flux,
    volume_integral = VolumeIntegralFluxDifferencing(volume_flux),
)

###############################################################################
# Create the TreeMesh for the domain [-3, 3]

coordinates_min = -3.0
coordinates_max = 3.0
mesh = TreeMesh(
    coordinates_min,
    coordinates_max,
    initial_refinement_level = 6,
    n_cells_max = 10_000,
    periodicity = true,
)

# create the semi discretization object
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver)

###############################################################################
# ODE solver

tspan = (0.0, 1.0)
ode = semidiscretize(semi, tspan)

###############################################################################
# Callbacks

summary_callback = SummaryCallback()

analysis_interval = 100
analysis_callback =
    AnalysisCallback(semi, interval = analysis_interval, save_analysis = false)

alive_callback = AliveCallback(analysis_interval = analysis_interval)

save_solution = SaveSolutionCallback(
    interval = 100,
    save_initial_solution = true,
    save_final_solution = true,
)

stepsize_callback = StepsizeCallback(cfl = 0.5)

callbacks = CallbackSet(
    summary_callback,
    analysis_callback,
    alive_callback,
    save_solution,
    stepsize_callback,
)

###############################################################################
# run the simulation

sol = solve(
    ode,
    CarpenterKennedy2N54(williamson_condition = false);
    dt = 1.0, # solve needs some value here but it will be overwritten by the stepsize_callback
    ode_default_options()...,
    callback = callbacks,
);
