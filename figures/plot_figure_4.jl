
using Trixi
using CairoMakie
using DataFrames
using DelimitedFiles

include("../src/main.jl")

## Run with WB configuration
trixi_include("../examples/elixir_shallowwater_moments_stone_throw.jl",
              equations = ShallowWaterMomentEquations1D(gravity = 9.812, H0 = 1.75,
                                                        n_moments = 2),
              surface_flux = (FluxPlusDissipation(flux_ec,
                                                  DissipationLaxFriedrichsEntropyVariables(Trixi.max_abs_speed)),
                              flux_nonconservative_ec))

# Get plot data
pd_wb = PlotData1D(sol)
# Get analysis data
mat, head = readdlm(joinpath("out", "analysis.dat"), header = true)
df_wb = DataFrame(mat, vec(head))

## Run with NWB configuration
trixi_include("../examples/elixir_shallowwater_moments_stone_throw.jl",
              equations = ShallowWaterMomentEquations1D(gravity = 9.812, H0 = 1.75,
                                                        n_moments = 2),
              surface_flux = (FluxPlusDissipation(flux_ec,
                                                  DissipationLocalLaxFriedrichs(Trixi.max_abs_speed)),
                              flux_nonconservative_ec))

# Get plot data
pd_nwb = PlotData1D(sol)
# Get analysis data
mat, head = readdlm(joinpath("out", "analysis.dat"), header = true)
df_nwb = DataFrame(mat, vec(head))

# Create plot
sizelblx = 20
sizelbly = 20

with_theme(theme_latexfonts()) do
    f = Figure(size = (800, 650))

    g_h = f[1, 1] = GridLayout()
    g_v = f[1, 2] = GridLayout()
    g_a1 = f[2, 1] = GridLayout()
    g_a2 = f[2, 2] = GridLayout()

    ax_h_wb = Axis(g_h[1, 1], xlabel = L"x", ylabel = L"h", xlabelsize = sizelblx,
                   ylabelsize = sizelbly)
    ax_v_wb = Axis(g_v[1, 1], xlabel = L"x", ylabel = L"u_m", xlabelsize = sizelblx,
                   ylabelsize = sizelbly)
    ax_a1_wb = Axis(g_a1[1, 1], xlabel = L"x", ylabel = L"\alpha_1", xlabelsize = sizelblx,
                    ylabelsize = sizelbly)
    ax_a2_wb = Axis(g_a2[1, 1], xlabel = L"x", ylabel = L"\alpha_2", xlabelsize = sizelblx,
                    ylabelsize = sizelbly)

    wb_color = Makie.wong_colors()[1]
    nwb_color = Makie.wong_colors()[2]

    label_wb = lines!(ax_h_wb, pd_wb.x, pd_wb.data[:, 1], color = wb_color, linewidth = 2)
    label_b = lines!(ax_h_wb, pd_wb.x, pd_wb.data[:, end], color = :black,
                     linestyle = :dash)
    lines!(ax_v_wb, pd_wb.x, pd_wb.data[:, 2], color = wb_color, linewidth = 2)
    lines!(ax_a1_wb, pd_wb.x, pd_wb.data[:, 3], color = wb_color, linewidth = 2)
    lines!(ax_a2_wb, pd_wb.x, pd_wb.data[:, 4], color = wb_color, linewidth = 2)

    label_nwb = lines!(ax_h_wb, pd_wb.x, pd_nwb.data[:, 1], color = nwb_color,
                       linewidth = 2)
    lines!(ax_v_wb, pd_wb.x, pd_nwb.data[:, 2], color = nwb_color, linewidth = 2)
    lines!(ax_a1_wb, pd_wb.x, pd_nwb.data[:, 3], color = nwb_color, linewidth = 2)
    lines!(ax_a2_wb, pd_wb.x, pd_nwb.data[:, 4], color = nwb_color, linewidth = 2)

    # Reset xlimits for all axes
    for ax in (ax_h_wb, ax_v_wb, ax_a1_wb, ax_a2_wb)
        Makie.xlims!(ax, (-4.0, 4.0))
        ax.xticklabelsize = 15
        ax.yticklabelsize = 15
    end

    legend = Legend(f[3, 1:2],
                    [label_wb, label_nwb, label_b],
                    ["WB", "NWB", L"b"],
                    labelsize = 20)

    legend.orientation = :horizontal

    save("Figure4.pdf", f)
end

## Figure for entropy over time and well-balanced error
with_theme(theme_latexfonts()) do
    f = Figure(size = (800, 350))
    ax_left = Axis(f[1, 1], xlabel = "t",
                   ylabel = L"\frac{1}{|\Omega|}\int_{\Omega} \mathbb{E}(t) \,\text{d}x", xlabelsize = 20, ylabelsize = 15,)
    ax_right = Axis(f[1, 2], xlabel = "t",
                    ylabel = L"\frac{1}{|\Omega|}\int_{\Omega} |h(t) - h(0)| \,\text{d}x", xlabelsize = 20,ylabelsize = 15,)

    # Plot entropy on left axis
    label_wb = lines!(ax_left, df_wb.time, df_wb.entropy, linewidth = 2)
    label_nwb = lines!(ax_left, df_nwb.time, df_nwb.entropy, linewidth = 2)

    # Plot well-balanced error on right axis
    lines!(ax_right, df_wb.time, df_wb.var"|H0-(h+b)|", linewidth = 2)
    lines!(ax_right, df_nwb.time, df_nwb.var"|H0-(h+b)|", linewidth = 2)

    legend = Legend(f[2, 1:2],
                    [label_wb, label_nwb],
                    ["WB", "NWB"], labelsize = 20)

    # Reset xlimits for all axes
    for ax in (ax_left, ax_right)
        Makie.xlims!(ax, (0.0, 8000.0))
        ax.xticklabelsize = 15
        ax.yticklabelsize = 15
    end

    legend.orientation = :horizontal

    save("Figure5.pdf", f)
end
