using DrWatson
@quickactivate "MOEAs"
using Metaheuristics
using HardTestProblems
using LinearAlgebra
using Plots; gr()

include("../src/MOEAs.jl")

function main()

    rw_problems = [
                   "pressure_vessel",
                   "vibrating_platform",
                   "two_bar_Truss_design_problems",
                   "weldan_beam_design"
                  ]
    max_runs = 31

    N = 100
    fcalls = 200_000
    iterations = 2fcalls*N
    for (fnum, p) in enumerate(rw_problems)
        best_res = nothing
        best_hv = 0
        best_seed = 0
        for seed in 1:max_runs
            @show fnum, seed
            options = Options(;iterations, seed, f_calls_limit=fcalls)
            f, conf =  get_RW_MOP_problem(p);
            bounds = Array([conf[:xmin] conf[:xmax]]')

            algorithm = MOEAs.CCMO_NSGAII(;N = N, options = options)
            res = optimize(f, bounds, algorithm)
            hv = PerformanceIndicators.hypervolume(res.population, conf[:nadir])
            if isnothing(best_res) || hv > best_hv
                best_res = res
                best_hv = hv
                @show best_hv
                best_seed = seed
            end
            
        end
        res = best_res

        results = Dict("F" => fvals(res.population),
                       "X" => positions(res.population),
                       "fmin" => ideal(res.population),
                       "fmax" => nadir(res.population),
                       "benchmark" => :application,
                       "fnum" => fnum,
                       "basealgorithm" => :CCMO,
                       "seed" => best_seed
                      )

        d = Dict("benchmark" => :application, "fnum" => fnum)
        wsave(datadir("approx_fronts/", savename(d, "bson")), results)
    end
end


main()
