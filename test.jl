using Plots; 
using DifferentialEquations
using Optimization
using OptimizationNLopt
using LinearAlgebra
using LaTeXStrings
specific_plot_kwargs = Dict{Symbol,Any}(
    # plotattr(:Plot)
    # plotattr(:Axis)
    :html_output_format => :png, # can try {:png,:svg} here
    :dpi => 300,                 # DPI=300 is minimum for good quality and print
    :grid => true,              # eliminate grid lines
    :gridstyle => :dash,         # if you actually want grid lines, choose from [:auto, :solid, :dash, :dot, :dashdot, :dashdotdot]
    :gridlinewidth => 0.5,       # set to 0 to make the grid and axes outlines disappear
    :framestyle => :box,         # add matlab-esque boundary box on plot
    :size => (450,450),          # typically (600,450) looks fine
    :topmargin => 0Plots.mm,
    :bottommargin => 0Plots.mm,
    :leftmargin => 0Plots.mm,
    :rightmargin => 5Plots.mm,
    :background_color_subplot => :white,
    :legend_background_color => nothing, # removes white background and adds translucence to legend, useful if you still want to see what's goin on behind a plot legend
    :legend_foreground_color => nothing, # removes outer black border of legend, just remove it alongside the white background of default legend to improve visibility
    :fontfamily => "computer modern", # DO NOT DEVIATE FROM THIS IF USING GR()
    :guidefontfamily => :match,       # match any fontfamily on title labels to axis labels
    :tickfontfamily => :match,       # match any fontfamily on title labels to ticks and tick labels
    :titlefontsize => 6, # for title fontsize
    :guidefontsize => 6, # for x and y labels
    :legendfontsize => 5, # for legend labels
    :tickfontsize => 6, # for axis ticks
    # :title => L"Sample text: $- \frac{{\hbar ^2}}{{2m}}\frac{{\partial ^2 \psi (x,t)}}{{\partial x^2 }} + U(x)\psi (x,t) = i\hbar \frac{{\partial \psi (x,t)}}{{\partial t}}$",
    # :xlabel => L"Sample text: $- \frac{{\hbar ^2}}{{2m}}\frac{{\partial ^2 \psi (x,t)}}{{\partial x^2 }} + U(x)\psi (x,t) = i\hbar \frac{{\partial \psi (x,t)}}{{\partial t}}$",
    # :ylabel => L"Sample text: $- \frac{{\hbar ^2}}{{2m}}\frac{{\partial ^2 \psi (x,t)}}{{\partial x^2 }} + U(x)\psi (x,t) = i\hbar \frac{{\partial \psi (x,t)}}{{\partial t}}$",

    # plotattr(:Series)
    # :label => L"Sample legend text: $a=\frac{{\hbar ^2}}{{2m}}$",
    :legend => :topleft,                # placement of legend, see `:(outer ?)(top/bottom ?)(right/left ?)`
    :linewidth => 1.0,
    :linestyle => :solid,               # see also [:auto, :solid, :dash, :dot, :dashdot, :dashdotdot]
    :linealpha => 1,
    # :linecolor => 1,
    # marker
    :markersize => 0.2,
    # :markershape => :star4,             # can choose from many, personal favorites = {:circle,:star4,:utriangle,:cross}
    :markeralpha => 1.0,
    :markercolor => nothing,
    # for that pesky marker outline on the markers
    :markerstrokewidth => 0.0,           # VERY IMPORTANT, make the marker outlines disappear -> 0
    :markerstrokestyle => :solid,        # the outline style on the marker
    :markerstrokecolor => nothing,       # can force this to be the same color as the marker, but still want outline
    :markerstrokealpha => 1,             # opacity of marker outline
)
gr(; specific_plot_kwargs...) # set gr backend

include("stateestimator.jl")
include("seoptimization.jl")




# Variables
tSpan = (0, 60)

Ca_0 = 20
Cb_0 = 30
Cc_0 = 0
k = 0.02


# True System
function reaction(du, u, p, t)
    du[1] = -k * u[1] * u[2]
    du[2] = -k * u[1] * u[2]
    du[3] = k * u[1] * u[2]
end


# State Estimator
L = [0.015068724372954935 1.0900011538659846e-22; 0.40464737008963436 0.00029748351739991287; 0.391685305678192 0.8731733818202637]
estX0 = [0,0,30]

trueProblem = ODEProblem(reaction, [Ca_0, Cb_0, Cc_0], tSpan)
estimatorProblem = ODEProblem(reaction, [0,0,0], tSpan)
function measurement(u)
    return [u[1] + u[2]; u[3]]
end

stateEstimator = StateEstimator(trueProblem, estimatorProblem, measurement, tSpan, Luenberger(L))
SE_PolyChaos(50, 2, 6).optimize(stateEstimator)
sol = solve(trueProblem, Tsit5())
estSol = solveEstimatedState(stateEstimator)


# Plot
function Save(dir) 
    savefig(string(@__DIR__, "/figs/$dir.png"))
end

graph = plot(sol, label=[L"$C_a$" L"$C_b$" L"$C_c$"], xlabel="Time (t)", ylabel="Concentration (mol)", title="Bilinear Reaction State Estimation")
plot!(graph, estSol, linestyle=:dash, label=[L"$C_a$ Estimation" L"$C_b$ Estimation" L"$C_c$ Estimation"], xlabel="Time (t)")
display(graph)

# Save("state-estimation-100x")