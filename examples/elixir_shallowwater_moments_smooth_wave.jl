
using OrdinaryDiffEqSSPRK, OrdinaryDiffEqLowStorageRK
using Trixi
using TrixiShallowWater

include("../src/main.jl")

# Semidiscretization of the shallow water moment equations

equations = ShallowWaterMomentEquations1D(gravity = 1.0, H0 = 0.0, n_moments = 2)

# Initial condition with a smooth wave wave in a periodic domain.
# See section 4.2 in the paper:
#   Julian Koellermeier and Marvin Rominger (2020)
#   "Analysis and Numerical Simulation of Hyperbolic Shallow Water Moment Equations"
#   [DOI: 10.4208/cicp.OA-2019-0065](https://doi.org/10.4208/cicp.OA-2019-0065)
function initial_condition_smooth_periodic_wave(x, t, equations::ShallowWaterMomentEquations1D)

    # Case JK 2020:
    H    = 1.0 + exp(3.0 * cos( π * (x[1] + 0.5)) )/exp(4.0)      
    v    = 0.25    
    a    = zeros(equations.n_moments)
    a[2] = -0.25
    
    b = 0.0

    return prim2cons(SVector(H, v, a..., b), equations)
end

initial_condition = initial_condition_smooth_periodic_wave

###############################################################################
# Get the DG approximation space

volume_flux  = (flux_ec, flux_nonconservative_ec)
surface_flux = (FluxPlusDissipation(flux_ec, DissipationLocalLaxFriedrichs(Trixi.max_abs_speed)), flux_nonconservative_ec)
#surface_flux = (FluxPlusDissipation(flux_ec, dissipation_pvm_force), flux_nonconservative_ec)

indicator_var(u, equations) = u[2]^3
basis = LobattoLegendreBasis(4)
indicator_sc = IndicatorHennemannGassner(equations, basis,
                                         alpha_max=0.5,
                                         alpha_min=0.001,
                                         alpha_smooth=true,
                                         variable=indicator_var)

volume_integral = VolumeIntegralShockCapturingHG(indicator_sc;
                                                 volume_flux_dg=volume_flux,
                                                 volume_flux_fv=surface_flux,)
#volume_integral = VolumeIntegralFluxDifferencing(volume_flux)

solver = DGSEM(basis, surface_flux, volume_integral)

###############################################################################
# Create the TreeMesh for the domain [-1, 1]

coordinates_min = -1.0
coordinates_max =  1.0

mesh = TreeMesh(coordinates_min,
                coordinates_max,
                initial_refinement_level = 7, # 2^refinement_level
                n_cells_max = 10_000,
                periodicity = true)           

# create the semi discretization object
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,  source_terms = source_term_bottom_friction)

###############################################################################
# ODE solver
tspan = (0.0, 2.0)
ode   = semidiscretize(semi, tspan)

###############################################################################
# Callbacks
summary_callback = SummaryCallback()

analysis_interval  =  10
analysis_callback  =  AnalysisCallback(semi, interval  =  analysis_interval, save_analysis  =  true,
                                       extra_analysis_integrals = (entropy,))
alive_callback     =  AliveCallback(analysis_interval  =  analysis_interval)
save_solution      =  SaveSolutionCallback(interval = 100, save_initial_solution = true, save_final_solution = true)
stepsize_callback  =  StepsizeCallback(cfl = 0.5)

callbacks = CallbackSet(summary_callback,  analysis_callback,  alive_callback,  save_solution,  stepsize_callback)
###############################################################################
# run the simulation

sol = solve(ode,
            CarpenterKennedy2N54(williamson_condition = false);
            dt = 1.0,              # solve needs some value here but it will be overwritten by the stepsize_callback
            ode_default_options()...,
            callback = callbacks,);
            
# CarpenterKennedy2N54: The five-stage, fourth order low-storage method            

###############################################################################
# Visualization
using CairoMakie

pd = PlotData1D(sol, reinterpolate=false)

with_theme(theme_latexfonts()) do
    f = Figure(size = (800, 600))

    g_h = f[1, 1] = GridLayout()
    g_v = f[1, 2] = GridLayout()
    g_a1 = f[2, 1] = GridLayout()
    g_a2 = f[2, 2] = GridLayout()
    
    ax_h = Axis(g_h[1,1], xlabel = L"x", ylabel = L"h")
    ax_v = Axis(g_v[1,1], xlabel = L"x", ylabel = L"u_m")
    ax_a1 = Axis(g_a1[1,1], xlabel = L"x", ylabel = L"\alpha_1")
    ax_a2 = Axis(g_a2[1,1], xlabel = L"x", ylabel = L"\alpha_2")

    lines!(ax_h, pd.x, pd.data[:,1])
    lines!(ax_v, pd.x, pd.data[:,2])
    lines!(ax_a1, pd.x, pd.data[:,3])
    lines!(ax_a2, pd.x, pd.data[:,4])

    # Reset xlimits for all axes
    for ax in (ax_h, ax_v, ax_a1, ax_a2)
        Makie.xlims!(ax, (-1, 1))
    end

    save("plot_smooth_wave.pdf", f)
end
