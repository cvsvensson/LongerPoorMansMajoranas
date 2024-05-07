
function plot_ss!(p, ss, N; kwargs...)
    ss_pos = ss.sol
    Plots.scatter!(p, [ss_pos[2]], [ss_pos[1]]; markersize=5, color=:red, kwargs...)
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

#https://discourse.julialang.org/t/tick-size-in-plots-jl/74793/3
function ticks_length!(; tl=0.02)
    p = Plots.current()
    xticks, yticks = Plots.xticks(p)[1][1], Plots.yticks(p)[1][1]
    xl, yl = Plots.xlims(p), Plots.ylims(p)
    x1, y1 = zero(yticks) .+ xl[1], zero(xticks) .+ yl[1]
    sz = p.attr[:size]
    r = sz[1] / sz[2]
    dx, dy = tl * (xl[2] - xl[1]), tl * r * (yl[2] - yl[1])
    Plots.plot!([xticks xticks]', [y1 y1 .+ dy]', c=:black, labels=false)
    Plots.plot!([x1 x1 .+ dx]', [yticks yticks]', c=:black, labels=false, xlims=xl, ylims=yl)
    return Plots.current()
end


function plot_f(data, f; clim_max=1, clims=missing, frame=:box, plot_ss=true, level_magnitude=0.01, ss_label="$(data["N"])-site sweet spot", kwargs...)
    x = data["x"]
    xlabel = data["labels"][1]
    y = data["y"]
    ylabel = data["labels"][2]
    levels = level_magnitude .* [-1, 0, 1]
    titles = [string(level_magnitude, "Δ"), "0"]
    colors = [:orange, :red, :orange]
    z = map(f, data["data"])

    xdiff = abs(first(x) - last(x)) / 2
    ydiff = abs(first(y) - last(y)) / 2
    xlims = (first(x), last(x)) #.+ xdiff / 40 .* (-1, 1)
    ylims = (first(y), last(y)) #.+ ydiff / 40 .* (-1, 1)

    clims = ismissing(clims) ? (0minimum(levels) - 1e-12, clim_max) : clims

    p = Plots.heatmap(x, y, z; clims, xlabel, ylabel, frame, xlims, ylims, kwargs...)
    ticks_length!(tl=0.015)
    zgap = map(x -> x.gap, data["data"])
    foreach((l, c) -> Plots.contour!(p, x, y, abs.(zgap); levels=[l], c), levels[[1, 3]], colors[[1, 3]])
    Plots.contour!(p, x, y, zgap; levels=[0.0], c=colors[2], clims)
    foreach((l, c) ->Plots.plot!(p, last(x) .* [1, 1] .+ 0.1, [0, 0]; leg_title="|δE|", c, label=l), titles, colors)
    #contour!(p, x, y, zgap; levels=range(-1, 1, 50))
    if plot_ss
        plot_ss!(p, data["ss"], data["N"]; label=ss_label)
    end
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
        xlims, ylims)
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