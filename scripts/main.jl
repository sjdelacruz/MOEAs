using DrWatson
@quickactivate "MOEAs"
using Metaheuristics
using HardTestProblems
using LinearAlgebra
using Plots; plotly()

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

        preferences = [ [0.1, 0.9], [0.5, 0.5], [0.8, 0.2] ]
        δ = fill(0.005, length(preferences))
        rw_problems = [
                       "pressure_vessel",
                       "vibrating_platform",
                       "two_bar_Truss_design_problems",
                       "weldan_beam_design"
                      ]
        # f, conf =  get_RW_MOP_problem(rw_problems[problem_idx]);
        # bounds = Array([conf[:xmin] conf[:xmax]]')
        f, bounds, _ = Metaheuristics.TestProblems.get_problem(:ZDT1)
        front = []
    end

    f, bounds, front, preferences, δ
end

function plot_res(res, preferences)
    p = plot(xlabel="f₁", ylabel="f₂", zlabel="f₃")
    fs = fvals(res)
    # normalize
    z_ideal = ideal(fs)'
    z_nad = nadir(fs)'
    #fs = (fs .- z_ideal) ./ (z_nad .- z_ideal)
    M = size(fs, 2)
    scatter!((fs[:,i] for i in 1:M)..., label="Result")

    for (i,w) in enumerate(preferences)
        t = range(0,1, length=50)
        # ww =  w .* (z_nad - z_ideal)# / norm(w)
        ww =  w # / norm(w)
        #line = z_ideal' .+ t.*ww'
        line = t.*ww'
        plot!((line[:,i] for i in 1:M)..., label="Preference $i", color=:blue, lw=2)
    end
    p
end

function get_ccmo_parms(d) 
    f, bounds, front, preferences, δ = get_problem(d[:fnum], d[:benchmark])
    options = Options(f_calls_limit=d[:fcalls],iterations=2d[:fcalls]*d[:N], seed=d[:seed], debug=false)

    algorithm = MOEAs.CCMO_NSGAII(;N = d[:N], options = options, preferences,δ)
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
    plot_res(res, algorithm.parameters.preferences)
end

main()

