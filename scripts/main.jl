using DrWatson
@quickactivate "MOEAs"
using Metaheuristics
using HardTestProblems
using LinearAlgebra
using Plots; gr()

include("../src/MOEAs.jl")

function run_algorithm(algorithm, f, bounds)
    res = Metaheuristics.optimize(f, bounds, algorithm)
    res
    
end


function get_problem(problem_idx, problem_type)

    # preferences defined as weights
    preferences = [[0.2,0.7, 0.1],
                   [0.25, 0.25, 0.5],
                   [0.7, 0.2, 0.1],
                  ]
    # the threshold (cosine distance)
    δ = fill(0.01, length(preferences))

    if problem_type == :constrianed        
        constrianed_problems =  [:C1_DTLZ1, :DC2_DTLZ2, :DC1_DTLZ3, :DC3_DTLZ4 ]
        problem = constrianed_problems[problem_idx]
        D = 7
        f, bounds, front = Metaheuristics.TestProblems.get_problem(problem)
        bounds = bounds[:, 1:D]
    elseif problem_type == :unconstrianed
        uncostrianed_problems = [:DTLZ1, :DTLZ2, :DTLZ3, :DTLZ4]
        problem = uncostrianed_problems[problem_idx]
        f, bounds, front = Metaheuristics.TestProblems.get_problem(problem)
    elseif problem_type ==:application

        weight_points = Vector{Float64}[
                       #[0.1, 0.9],
                       #[0.5, 0.5],
                       # [0.8, 0.2]
                      ]
        ref_points = Vector{Float64}[
                      [-0.032, 260]
                     ]
        δ_w = fill(0.1, length(weight_points))
        δ_r = fill(0.01, length(ref_points))
        rw_problems = [
                       "pressure_vessel",
                       "vibrating_platform",
                       "two_bar_Truss_design_problems",
                       "weldan_beam_design"
                      ]
        f, conf =  get_RW_MOP_problem(rw_problems[problem_idx]);
        bounds = Array([conf[:xmin] conf[:xmax]]')
        # f, bounds, _ = Metaheuristics.TestProblems.get_problem(:ZDT3)
        front = []
    end

    f, bounds, front, weight_points, δ_w, ref_points, δ_r
end

function plot_res(archive, population, weight_points, ref_points)
    M = length(Metaheuristics.fval(population[1]))
    p = plot(xlabel="f₁", ylabel="f₂", zlabel="f₃", dpi=200)

    # plot ref directions
    z_ideal = ideal(population)
    z_nad = nadir(population)
    for (i,w) in enumerate(weight_points)
        t = range(0,1, length=50)
        # scale to axis
        ww =  w .* (z_nad - z_ideal)
        line = z_ideal' .+ t.*ww'

        plot!((line[:,i] for i in 1:M)..., label="", color=:gray, lw=2)
    end


    for r in ref_points
        scatter!(r[1:1], r[2:2], markercolor=:red, label="")
    end

    # population
    fs = fvals(population)
    scatter!((fs[:,i] for i in 1:M)..., label="Approx. Front", markercolor=:lightgray, markerstrokewidth=0)

    fs = fvals(archive)
    n = length(archive)
    scatter!((fs[:,i] for i in 1:M)..., label="Prefered solutions ($n)", markercolor=:black)
    savefig(p, "fig.png")


    p
end

function get_ccmo_parms(d) 

    f, bounds, front, weight_points, δ_w, ref_points, δ_r = get_problem(d[:fnum], d[:benchmark])
    options = Options(f_calls_limit=d[:fcalls],iterations=2d[:fcalls]*d[:N], seed=d[:seed], debug=false)

    algorithm_ = MOEAs.CCMO_NSGAII(;N = d[:N], options = options)
    algorithm = MOEAs.ROIs(algorithm_; weight_points, δ_w, ref_points, δ_r)
    algorithm, f, bounds
end


function main()
    seed = 1
    fnum = 2
    benchmark = :application

    d = Dict(
             :fnum => fnum,
             :benchmark => benchmark,
             :fcalls => 100_000,
             :N => 100,
             :seed => seed
            )

    algorithm, f, bounds = get_ccmo_parms(d)
    res = run_algorithm(algorithm, f, bounds)
    display(res)
    display(algorithm.parameters.archive)
    plot_res(algorithm.parameters.archive, res.population, algorithm.parameters.weight_points,algorithm.parameters.ref_points)
end

main()

