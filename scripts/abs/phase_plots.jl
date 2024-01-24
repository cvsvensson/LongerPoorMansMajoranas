function plot_LD(data)
    x = data["x"]
    xlabel = data["labels"][1]
    y = data["y"]
    ylabel = data["labels"][2]
    ss = data["ss"]
    ss_pos = ss.sol
    levels = 0.01 .* [-1, 1, 0]
    z = map(LDf, data["data"])
    p = heatmap(x, y, z; c=cgrad(:viridis, rev=true), clims=(minimum(levels) - 1e-12, 1), xlabel, ylabel, colorbar_title="LD")
    colors = [:white, :white, :red]
    zgap = map(x -> x.gap, data["data"])
    foreach((l, c) -> contour!(p, x, y, zgap; levels=[l], c), levels, colors)
    scatter!(p, [ss_pos[2]], [ss_pos[1]], markersize=5, color=:red, label="$(data["N"])-site sweet spot")
    p
end
function plot_gap(data)
    x = data["x"]
    xlabel = data["labels"][1]
    y = data["y"]
    ylabel = data["labels"][2]
    ss = data["ss"]
    ss_pos = ss.sol
    levels = 0.01 .* [-1, 1, 0]
    z = map(x -> tanh(x.gap), data["data"])
    p = heatmap(x, y, z; c=:redsblues, clims=(-1, 1), xlabel, ylabel)
    colors = [:gray, :gray, :black]
    zgap = map(x -> x.gap, data["data"])
    foreach((l, c) -> contour!(p, x, y, zgap; levels=[l], c), levels, colors)
    scatter!(p, [ss_pos[2]], [ss_pos[1]], markersize=5, color=:red, label="$(data["N"])-site sweet spot")
    p
end
function plot_MPU(data)
    x = data["x"]
    xlabel = data["labels"][1]
    y = data["y"]
    ylabel = data["labels"][2]
    ss = data["ss"]
    ss_pos = ss.sol
    levels = 0.01 .* [-1, 1, 0]
    z = map(MPU, data["data"])
    p = heatmap(x, y, z; c=cgrad(:viridis, rev=true), clims=(minimum(levels) - 1e-12, 1), xlabel, ylabel, colorbar_title="MPU")
    colors = [:white, :white, :red]
    zgap = map(x -> x.gap, data["data"])
    foreach((l, c) -> contour!(p, x, y, zgap; levels=[l], c), levels, colors)
    scatter!(p, [ss_pos[2]], [ss_pos[1]], markersize=5, color=:red, label="$(data["N"])-site sweet spot")
    p
end
function plot_MP(data)
    x = data["x"]
    xlabel = data["labels"][1]
    y = data["y"]
    ylabel = data["labels"][2]
    ss = data["ss"]
    ss_pos = ss.sol
    levels = 0.01 .* [-1, 1, 0]
    z = map(MP, data["data"])
    p = heatmap(x, y, z; c=cgrad(:viridis, rev=true), clims=(minimum(levels) - 1e-12, 1), xlabel, ylabel, colorbar_title="MP")
    colors = [:white, :white, :red]
    zgap = map(x -> x.gap, data["data"])
    foreach((l, c) -> contour!(p, x, y, zgap; levels=[l], c), levels, colors)
    scatter!(p, [ss_pos[2]], [ss_pos[1]], markersize=5, color=:red, label="$(data["N"])-site sweet spot")
    p
end