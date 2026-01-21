
using OrdinaryDiffEqSSPRK, OrdinaryDiffEqLowStorageRK
using Trixi
using TrixiShallowWater

include("../src/main.jl")

# Semidiscretization of the shallow water moment equations
equations = ShallowWaterLinearizedMomentEquations1D(gravity = 9.812, H0 = 1.75, n_moments = 3)

function initial_condition_stone_throw(x, t, equations::ShallowWaterLinearizedMomentEquations1D)
    H    = 1.75
    
    # Set discontinuous velocity
    v    = 0.0   
    a    = zeros(equations.n_moments)
    if x[1] >= -0.75 && x[1] <= 0.0
        v = -1.0
    elseif x[1] > 0.0 && x[1] <= 0.75
        v = 1.0
    end

    # Set smooth bottom topography
    b = 1 / exp(0.5 * x[1]^2)

    h = H - b

    return SVector(h, h*v, (h*a)..., b)
end

initial_condition = initial_condition_stone_throw

###############################################################################
# Get the DG approximation space

volume_flux  = (flux_ec, flux_nonconservative_ec)
surface_flux = (FluxPlusDissipation(flux_ec, DissipationLaxFriedrichsEntropyVariables(Trixi.max_abs_speed)), flux_nonconservative_ec)

indicator_var(u, equations) = u[2]^3
basis = LobattoLegendreBasis(3)
indicator_sc = IndicatorHennemannGassner(equations, basis,
                                         alpha_max=0.5,
                                         alpha_min=0.001,
                                         alpha_smooth=true,
                                         variable=indicator_var)

volume_integral = VolumeIntegralShockCapturingHG(indicator_sc;
                                                 volume_flux_dg=volume_flux,
                                                 volume_flux_fv=surface_flux,)

solver = DGSEM(basis, surface_flux, volume_integral)

###############################################################################
# Create the TreeMesh for the domain [-1, 1]

coordinates_min = -3.0
coordinates_max =  3.0

mesh = TreeMesh(coordinates_min,
                coordinates_max,
                initial_refinement_level = 5, # 2^refinement_level
                n_cells_max = 10_000,
                periodicity = true)           

# create the semi discretization object
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver)
#,  source_terms = source_term_bottom_friction

###############################################################################
# ODE solver
tspan = (0.0, 100.0)
ode   = semidiscretize(semi, tspan)

###############################################################################
# Callbacks
summary_callback = SummaryCallback()

# @inline function dwdP(u, equations)
#     w = cons2entropy(u, equations)
#     P = source_term_bottom_friction(u, zero(eltype(u)), zero(eltype(u)), equations)

#     return w' * P
# end

analysis_interval  =  100000
analysis_callback  =  AnalysisCallback(semi, interval  =  analysis_interval, save_analysis  =  true,
                                       extra_analysis_integrals = (entropy, lake_at_rest_error),)
alive_callback     =  AliveCallback(analysis_interval  =  analysis_interval)
save_solution      =  SaveSolutionCallback(dt = 100.0, save_initial_solution = true, save_final_solution = true)
stepsize_callback  =  StepsizeCallback(cfl = 0.9)

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
# using CairoMakie

# pd = PlotData1D(sol, reinterpolate=false)

# with_theme(theme_latexfonts()) do
#     f = Figure(size = (800, 600))

#     g_h = f[1, 1] = GridLayout()
#     g_v = f[1, 2] = GridLayout()
#     g_a1 = f[2, 1] = GridLayout()
#     g_a2 = f[2, 2] = GridLayout()
#     g_a3 = f[3, 1] = GridLayout()
    
#     ax_h = Axis(g_h[1,1], xlabel = L"x", ylabel = L"h + b")
#     ax_v = Axis(g_v[1,1], xlabel = L"x", ylabel = L"u_m")
#     ax_a1 = Axis(g_a1[1,1], xlabel = L"x", ylabel = L"\alpha_1")
#     ax_a2 = Axis(g_a2[1,1], xlabel = L"x", ylabel = L"\alpha_2")
#     ax_a3 = Axis(g_a3[1,1], xlabel = L"x", ylabel = L"\alpha_3")

#     lines!(ax_h, pd.x, pd.data[:,1])
#     lines!(ax_h, pd.x, pd.data[:,end])
#     lines!(ax_v, pd.x, pd.data[:,2])
#     lines!(ax_a1, pd.x, pd.data[:,3])
#     lines!(ax_a2, pd.x, pd.data[:,4])
#     #lines!(ax_a3, pd.x, pd.data[:,5])

#     # Reset xlimits for all axes
#     for ax in (ax_h, ax_v, ax_a1, ax_a2, ax_a3)
#         Makie.xlims!(ax, (-3.0, 3.0))
#     end

#     save("plot_stone_throw.pdf", f)
# end
