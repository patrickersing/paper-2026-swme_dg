
using Trixi
using CairoMakie
using DelimitedFiles
using DataFrames
using LaTeXStrings

include("../src/main.jl")

pd_list = []
ref_levels = [7, 8, 9]

# Run examples for different polynomial degrees and save plot data
for ref_level in ref_levels
    trixi_include("../examples/elixir_shallowwater_moments_smooth_wave.jl",
                  initial_refinement_level = ref_level)
    
    pd = PlotData1D(sol)
    push!(pd_list, pd)
end

# Read reference file
pd_ref = readdlm("../figures/soln10_Euler_HLL2500.dat", Float64, header = false)

# Create plot
colors = Makie.wong_colors()
sizelblx = 25
sizelbly = 25

with_theme(theme_latexfonts()) do
    f = Figure(size = (800, 650))

    g_h = f[1, 1] = GridLayout()
    g_v = f[1, 2] = GridLayout()
    g_a1 = f[2, 1] = GridLayout()
    g_a2 = f[2, 2] = GridLayout()
    
    ax_h = Axis(g_h[1,1], xlabel = L"x", ylabel = L"h", xlabelsize = sizelblx, ylabelsize = sizelbly)
    ax_v = Axis(g_v[1,1], xlabel = L"x", ylabel = L"u_m", xlabelsize = sizelblx, ylabelsize = sizelbly)
    ax_a1 = Axis(g_a1[1,1], xlabel = L"x", ylabel = L"\alpha_1", xlabelsize = sizelblx, ylabelsize = sizelbly)
    ax_a2 = Axis(g_a2[1,1], xlabel = L"x", ylabel = L"\alpha_2", xlabelsize = sizelblx, ylabelsize = sizelbly)

    # Plot data for polynomial degrees
    for (i, ref_level) in enumerate(ref_levels)
        pd = pd_list[i]
        lines!(ax_h, pd.x, pd.data[:,1], color = colors[i], linewidth = 2)
        lines!(ax_v, pd.x, pd.data[:,2], color = colors[i], linewidth = 2)
        lines!(ax_a1, pd.x, pd.data[:,3], color = colors[i], linewidth = 2)
        lines!(ax_a2, pd.x, pd.data[:,4], color = colors[i], linewidth = 2)
    end

    # Add reference
    lines!(ax_h,  pd_ref[:,1], pd_ref[:,2], label = "Ref", color = :black, linestyle = :dash, linewidth = 3)    
    lines!(ax_v,  pd_ref[:,1], pd_ref[:,3], label = "Ref", color = :black, linestyle = :dash, linewidth = 3)
    lines!(ax_a1, pd_ref[:,1], pd_ref[:,4], label = "Ref", color = :black, linestyle = :dash, linewidth = 3)
    lines!(ax_a2, pd_ref[:,1], pd_ref[:,5], label = "Ref", color = :black, linestyle = :dash, linewidth = 3)

    # Reset xlimits for all axes
    for ax in (ax_h, ax_v, ax_a1, ax_a2)
        Makie.xlims!(ax, (-1.0, 1.0))
        ax.xticklabelsize = 18
        ax.yticklabelsize = 18
    end

    # Add legend
    legend = Legend(f[3,1:2], 
        [[LineElement(linestyle = :solid, color = colors[i], linewidth = 2) for i in 1:length(ref_levels)]...,
         LineElement(linestyle = :dash, color = :black, linewidth = 3)],
         [[L"K=%$(2^ref_level)" for ref_level in ref_levels]..., "Ref"],
         labelsize = 18)

    legend.orientation = :horizontal

    save("Figure2.pdf", f)
end