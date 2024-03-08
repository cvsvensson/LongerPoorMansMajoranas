
default(; c=cgrad(:viridis, rev=true, scale=x -> exp(5abs(x))),
    fontfamily="Computer Modern",
    colorbar_titlefontrotation=-90
)
##
Ns = [2, 3, 4, 5, 20, 40]
Fdatas = Dict(zip(Ns, load_full_data.(Ns)))
Fdatas_ss = Dict(zip(Ns, load_full_data.(Ns[1:4], "ss")))
Fdatas_ss2 = Dict(zip(Ns, load_full_data.(Ns[1:4], "ss")))
foreach(N -> (Fdatas[N]["ss"] = Fdatas_ss[N]["ss"]), Ns[1:4])
foreach(N -> (Fdatas[N]["ss"] = Fdatas_ss2[N]["ss"]), Ns[1:4])
##
phase_plots = [plot_f(Fdatas[N], MP; leg=:bottomleft, plot_ss=true, title="$N sites", colorbar_title="1-MP", ss_label=false) for N in Ns]
##
display.(phase_plots)
##
foreach((p, N) -> save(plotsdir(string("phase_diagram_MP_$N", ".png")), p, px_per_unit=4), phase_plots, Ns)


##
plot_f(F2data, MP, leg=:bottomleft, plot_ss=false, title="2 sites", colorbar_title="1-MP")
plot_f(F3data, MP, leg=:bottomleft, plot_ss=false, title="3 sites", colorbar_title="1-MP")
plot_f(F4data, MP, leg=:bottomleft, plot_ss=false, title="4 sites", colorbar_title="1-MP")
plot_f(F5data, MP, leg=:bottomleft, plot_ss=false, title="5 sites", colorbar_title="1-MP")
plot_f(F20data, MP, leg=:bottomleft, plot_ss=false, title="20 sites", colorbar_title="1-MP")
plot_f(F40data, MP, leg=:bottomleft, plot_ss=false, title="40 sites", colorbar_title="1-MP")

##
p2 = plot_f(F2data, MP, leg=:bottomleft, colorbar=false, plot_ss=false, title="2 sites", colorbar_title="1-MP", aspect_ratio=0.15)
p40 = plot_f(F40data, MP, leg=false, ylabel="", yshowaxis=false, colorbar=true, plot_ss=false, title="40 sites", right_margin=0Plots.mm, left_margin=0Plots.mm, colorbar_title="1-MP", aspect_ratio=0.16)
plot_ss!(p2, F2data["ss"], F2data["N"])
plot_ss!(p40, F2data["ss"], F2data["N"])
plot(p2, p40, size=(1000, 400), thickness_scaling=1.4)

##
plot_f(F4data, MP; c=cgrad(:viridis, rev=true, scale=x -> exp(4x)), plot_ss=false, fontfamily="Computer Modern", title="4", colorbar_titlefontrotation=-90, leg=:bottomleft, colorbar_title="1-MP")
plot_f(F5data, MP; c=cgrad(:viridis, rev=true, scale=x -> exp(4x)), plot_ss=false, fontfamily="Computer Modern", title="5", colorbar_titlefontrotation=-90, leg=:bottomleft, colorbar_title="1-MP")
plot_f(F20data, MP; c=cgrad(:viridis, rev=true, scale=x -> exp(4x)), plot_ss=false, fontfamily="Computer Modern", title="20", colorbar_titlefontrotation=-90, leg=:bottomleft, colorbar_title="1-MP")
plot_f(F40data, MP; c=cgrad(:viridis, rev=true, scale=x -> exp(6x)), plot_ss=false, fontfamily="Computer Modern", title="40", colorbar_titlefontrotation=-90, leg=:bottomleft, colorbar_title="1-MP")
