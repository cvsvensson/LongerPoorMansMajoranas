
function plot_ss!(p, ss, N)
    ss_pos = ss.sol
    scatter!(p, [ss_pos[2]], [ss_pos[1]], markersize=5, color=:red, label="$(N)-site sweet spot")
end
plot_ss!(p, ss::Union{Nothing,Missing}, N) = nothing
function plot_LD(data; clims=missing)
    x = data["x"]
    xlabel = data["labels"][1]
    y = data["y"]
    ylabel = data["labels"][2]
    levels = 0.01 .* [-1, 1, 0]
    z = map(LDf, data["data"])

    clims = ismissin(clims) ? (minimum(levels) - 1e-12, 1) : clims
    p = heatmap(x, y, z; c=cgrad(:viridis, rev=true), clims, xlabel, ylabel, colorbar_title="LD", frame=:box)
    colors = [:white, :white, :red]
    zgap = map(x -> x.gap, data["data"])
    plot_ss!(p, data["ss"], data["N"])
    foreach((l, c) -> contour!(p, x, y, zgap; levels=[l], c), levels, colors)
    p
end
function plot_gap(data)
    x = data["x"]
    xlabel = data["labels"][1]
    y = data["y"]
    ylabel = data["labels"][2]
    levels = 0.01 .* [-1, 1, 0]
    z = map(x -> tanh(x.gap), data["data"])
    p = heatmap(x, y, z; c=:redsblues, clims=(-1, 1), xlabel, ylabel, frame=:box)
    colors = [:gray, :gray, :black]
    zgap = map(x -> x.gap, data["data"])
    foreach((l, c) -> contour!(p, x, y, zgap; levels=[l], c), levels, colors)
    plot_ss!(p, data["ss"], data["N"])
    p
end
function plot_f(data, f; colorbar_title="", clims=missing)
    x = data["x"]
    xlabel = data["labels"][1]
    y = data["y"]
    ylabel = data["labels"][2]
    levels = 0.01 .* [-1, 1, 0]
    z = map(f, data["data"])

    xlims = (first(x), last(x))
    ylims = (first(y), last(y))

    clims = ismissing(clims) ? (minimum(levels) - 1e-12, 1) : clims

    p = heatmap(x, y, z; c=cgrad(:viridis, rev=true), clims, xlabel, ylabel, colorbar_title, frame=:box, colorbar_titlefontrotation=-90,
        xlims, ylims, thickness_scaling=1.3)
    colors = [:white, :white, :red]
    zgap = map(x -> x.gap, data["data"])
    foreach((l, c) -> contour!(p, x, y, zgap; levels=[l], c), levels, colors)
    plot_ss!(p, data["ss"], data["N"])
    p
end
function plot_MPU(data)
    x = data["x"]
    xlabel = data["labels"][1]
    y = data["y"]
    ylabel = data["labels"][2]
    levels = 0.01 .* [-1, 1, 0]
    z = map(MPU, data["data"])

    xlims = (first(x), last(x))
    ylims = (first(y), last(y))

    p = heatmap(x, y, z; c=cgrad(:viridis, rev=true), clims=(minimum(levels) - 1e-12, 1), xlabel, ylabel, colorbar_title="MPU", frame=:box, colorbar_titlefontrotation=-90,
        xlims, ylims, thickness_scaling=1.3)
    colors = [:white, :white, :red]
    zgap = map(x -> x.gap, data["data"])
    foreach((l, c) -> contour!(p, x, y, zgap; levels=[l], c), levels, colors)
    plot_ss!(p, data["ss"], data["N"])
    p
end
function plot_MP(data)
    x = data["x"]
    xlabel = data["labels"][1]
    y = data["y"]
    ylabel = data["labels"][2]
    levels = 0.01 .* [-1, 1, 0]
    z = map(MP, data["data"])
    p = heatmap(x, y, z; c=cgrad(:viridis, rev=true), clims=(minimum(levels) - 1e-12, 1), xlabel, ylabel, colorbar_title="MP")
    colors = [:white, :white, :red]
    zgap = map(x -> x.gap, data["data"])
    foreach((l, c) -> contour!(p, x, y, zgap; levels=[l], c), levels, colors)
    plot_ss!(p, data["ss"], data["N"])
    p
end