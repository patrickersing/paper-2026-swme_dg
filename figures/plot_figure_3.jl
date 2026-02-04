
using Trixi
using CairoMakie
using DelimitedFiles
using DataFrames
using LaTeXStrings

include("../src/main.jl")

df_ns = []
df_mm = []
nu_vals = [0.1, 0.2, 0.5, 1.0]

# Run examples with newtonian slip friction
for nu_val in nu_vals
    trixi_include("../examples/elixir_shallowwater_moments_smooth_wave.jl",
                  equations = ShallowWaterLinearizedMomentEquations1D(gravity = 9.81,
                                                                      H0 = 0.0,
                                                                      n_moments = 2,
                                                                      nu = nu_val),
                  source_terms = source_term_bottom_friction,
                  extra_analysis_integrals = (dwdP_Ns,))

    mat, head = readdlm(joinpath("out", "analysis.dat"), header = true)
    push!(df_ns, DataFrame(mat, vec(head)))
end

# Run examples with manning friction
for nu_val in nu_vals
    trixi_include("../examples/elixir_shallowwater_moments_smooth_wave.jl",
                  equations = ShallowWaterLinearizedMomentEquations1D(gravity = 9.81,
                                                                      H0 = 0.0,
                                                                      n_moments = 2,
                                                                      nu = nu_val),
                  source_terms = source_term_manning_friction,
                  extra_analysis_integrals = (dwdP_MM,))

    mat, head = readdlm(joinpath("out", "analysis.dat"), header = true)
    push!(df_mm, DataFrame(mat, vec(head)))
end

# Create plot
colors = Makie.wong_colors()
sizelblx = 20
sizelbly = 20

with_theme(theme_latexfonts()) do
    f = Figure(size = (800, 400))

    g_ns = f[1, 1] = GridLayout()
    g_mm = f[1, 2] = GridLayout()

    ax_ns = Axis(g_ns[1, 1], xlabel = L"t", ylabel = L"\mathcal{D}_{\textrm{Ns}}",
                 xlabelsize = sizelblx, ylabelsize = sizelbly)
    ax_mm = Axis(g_mm[1, 1], xlabel = L"t", ylabel = L"\mathcal{D}_{\textrm{NM}}",
                 xlabelsize = sizelblx, ylabelsize = sizelbly)

    # Plot data for polynomial degrees
    for (i, nu_val) in enumerate(nu_vals)
        lines!(ax_ns, df_ns[i].time, df_ns[i].var"dwdP_Ns", color = colors[i],
               linewidth = 2)
        lines!(ax_mm, df_mm[i].time, df_mm[i].var"dwdP_MM", color = colors[i],
               linewidth = 2)
    end

    # Plot reference line at zero
    lines!(ax_ns, df_ns[1].time, 0.0 * df_ns[1].var"dwdP_Ns", color = :black,
           linestyle = :dash, linewidth = 1)
    lines!(ax_mm, df_mm[1].time, 0.0 * df_mm[1].var"dwdP_MM", color = :black,
           linestyle = :dash, linewidth = 1)

    # Reset xlimits for all axes
    for ax in (ax_ns, ax_mm)
        Makie.xlims!(ax, (-0.01, 1.0))
        ax.xticklabelsize = 15
        ax.yticklabelsize = 15
    end

    # Add legend
    legend = Legend(f[2, 1:2],
                    [[LineElement(linestyle = :solid, color = colors[i], linewidth = 2)
                      for i in 1:length(nu_vals)]...],
                    [[L"ν=%$(nu_val)" for nu_val in nu_vals]...],
                    labelsize = 20)

    legend.orientation = :horizontal

    save("Figure3.pdf", f)
end
