using DrWatson
@quickactivate "MOEAs"
using Metaheuristics
using HardTestProblems
using LinearAlgebra
using BSON
using Plots; gr()

include("../src/MOEAs.jl")

function run_algorithm(algorithm, f, bounds)
    res = Metaheuristics.optimize(f, bounds, algorithm)
    res
    
end


function get_problem(problem_idx, problem_type, m = 2)

    # preferences defined as weights
    weight_points = [[0.2,0.7, 0.1],
                   [0.25, 0.25, 0.5],
                   [0.7, 0.2, 0.1],
                  ]
    # the threshold (cosine distance)
    δ_w = fill(0.2, length(weight_points))

    if m == 2
        weight_points = Vector{Float64}[
                                        [0.1, 0.9],
                                        [0.5, 0.5],
                                        [0.8, 0.2]
                                       ]
        δ_w = fill(0.2, length(weight_points))
    elseif m == 3
        # preferences defined as weights
        weight_points = [[0.2,0.7, 0.1],
                         [0.25, 0.25, 0.5],
                         [0.7, 0.2, 0.1],
                        ]
        # the threshold (cosine distance)
        δ_w = fill(0.2, length(weight_points))
    else
        error("only for m = 2,3")
    end
    

    if problem_type == :constrianed        
        constrianed_problems =  [:C1_DTLZ1, :C2_DTLZ2, :C1_DTLZ3, :C3_DTLZ4 ]
        problem = constrianed_problems[problem_idx]
        f, bounds, front = Metaheuristics.TestProblems.get_problem(problem)

        if m == 3
            D = 7
            bounds = bounds[:, 1:D]
        end
        expr = Meta.parse("Metaheuristics.TestProblems."* string(problem) )
        getp = eval(expr)
        ff, bounds, front = getp(m, gen_ref_dirs(m,100)) 
        f = x -> ff(x, m)
    elseif problem_type == :unconstrianed
        uncostrianed_problems = [:DTLZ1, :DTLZ2, :DTLZ3, :DTLZ4]
        problem = uncostrianed_problems[problem_idx]
        expr = Meta.parse("Metaheuristics.TestProblems."* string(problem) )
        getp = eval(expr)
        ff, bounds, front = getp(m, gen_ref_dirs(m,100)) 
        f = x -> ff(x, m)
    elseif problem_type ==:application

        rw_problems = [
                       "pressure_vessel",
                       "vibrating_platform",
                       "two_bar_Truss_design_problems",
                       "weldan_beam_design"
                      ]
        f, conf =  get_RW_MOP_problem(rw_problems[problem_idx]);
        bounds = Array([conf[:xmin] conf[:xmax]]')
        data = BSON.load(joinpath(scriptsdir(), "approx_fronts", "benchmark=application_fnum=$problem_idx.bson"))
        fs = data["F"]
        front = [ Metaheuristics.create_child(zeros(0), (fs[i,:],[0.0],[0.0])) for i in 1:size(fs, 1) ]
    end

    f, bounds, front, weight_points, δ_w
end

function plot_res(archive, res, d)
    population = res.population
    f, bounds, pf, weight_points, δ_w = get_problem(d["fnum"], d["benchmark"])

    M = length(Metaheuristics.fval(population[1]))
    p = plot(xlabel="f₁", ylabel="f₂", zlabel="f₃", dpi=200)

    if isempty(pf)
        z_ideal = ideal(population)
        z_nad = nadir(population)
    else
        z_ideal = ideal(pf)
        z_nad = nadir(pf)
    end
 
    # plot ref directions
    for (i,w) in enumerate(weight_points)
        t = range(0,1, length=50)
        if M <= 3
            w = w ./ norm(w)
        end
        
        # scale to axis
        ww =  w .* (z_nad - z_ideal)
        line = z_ideal' .+ t.*ww'


        if M == 2
            δ = δ_w[i]/2 # .* t .* norm(z_nad - z_ideal)
            #ww = (w .- δ) .* (z_nad - z_ideal)
            θ = max(atan(w[2] / w[1]) - π*δ/2,  0)
            ww = [cos(θ), sin(θ)] .* (z_nad - z_ideal)
            ln = z_ideal' .+ t.* ww'
            #plot!(ln[:,1], ln[:,2] .+ δ,  fill = (line[:,2] .- δ, :lightgray, 0.5), primary = false,linecolor = nothing,)
            #
            plot!((ln[:,i] for i in 1:M)..., label="", color=:gray, ls = :dot)
            θ = min(atan(w[2] / w[1]) + π*δ/2,  π/2)
            ww = [cos(θ), sin(θ)] .* (z_nad - z_ideal)
            ln = z_ideal' .+ t.* ww'
            plot!((ln[:,i] for i in 1:M)..., label="", color=:gray, ls = :dot)
        end
        plot!((line[:,i] for i in 1:M)..., label="", color=:gray, lw=2, arrow=true)
        
    end



    #= population
    fs = fvals(population)
    scatter!((fs[:,i] for i in 1:M)..., label="Approx. Front", markercolor=:lightgray, markerstrokewidth=0)
    =#


    if !isempty(pf)
        fs = fvals(pf)
        if size(fs, 2) == 2
            mask = sortperm(fs[:,1])
            plot!(fs[mask,1], fs[mask,2], label="Pareto Front", color = :darkgray)
        else
            scatter!((fs[:,i] for i in 1:M)...,
                     label="Pareto Front",
                     markercolor=:lightgray,
                     markerstrokewidth=0)
        end
        
    end

    fs = fvals(archive)
    scatter!((fs[:,i] for i in 1:M)..., label="Not in ROI", markercolor=:white)

    # identify those solutions in ROI
    mask = get_roi_perm(archive, z_ideal, z_nad, weight_points, δ_w)
    if  !isempty(mask)
        archive_ok = archive[mask]
        fs = fvals(archive_ok)
        scatter!((fs[:,i] for i in 1:M)..., label="In ROI", markercolor=:black)
    end

    p
end

function get_roi_perm(archive, fmin, fmax, w, δ)
    δ = 1 .- cos.(π/2*(δ/2)) # re-scale
    mask = MOEAs.compute_rio_vio(archive, fmin, fmax, w, δ) .<= 0
    return findall(mask)
end


function get_parms(d) 

    N = d["N"]
    f, bounds, front, weight_points, δ_w = get_problem(d["fnum"], d["benchmark"])
    options = Options(f_calls_limit=d["fcalls"],iterations=2d["fcalls"]*N, seed=d["seed"], debug=false)

    if d["basealgorithm"] == :CCMO
        base_algorithm = MOEAs.CCMO_NSGAII(;N = N, options = options)
    elseif d["basealgorithm"] == :SPEA2
        base_algorithm = SPEA2(;N , options = options)
    elseif d["basealgorithm"] == :SMS_EMOA
        base_algorithm = SMS_EMOA(;N , options = options)
    elseif d["basealgorithm"] == :NSGA2
        base_algorithm = NSGA2(;N, options = options)
        #=
        elseif d["basealgorithm"] == :MOEADDE
        nobjectives = d["benchmark"] in [:application] ? 2 : 3
        npartitions = d["benchmark"] in [:application] ? 100 : 12
        weights = gen_ref_dirs(nobjectives, npartitions)
        base_algorithm = MOEAD_DE(weights; options = options)
        =#
    else
        alg = d["basealgorithm"]
        error("Base algorithm $alg is not supported.")
    end
    
    algorithm = MOEAs.ROIs(base_algorithm; weight_points, δ_w)
    algorithm, f, bounds
end

function compute_roi_hv(archive, d)
    if isempty(archive)
        return 0.0
    end

    f, bounds, front, weight_points, δ_w = get_problem(d["fnum"], d["benchmark"])

    fmin = ideal(front)
    fmax = nadir(front)

    # keep only the ROI
    mask = get_roi_perm(archive, fmin, fmax, weight_points, δ_w)
    if isempty(mask)
        return 0.0
    end

    return PerformanceIndicators.hypervolume(archive[mask], fmax)
    
end


function main()
    seed = collect(1:31)
    fnum = collect(1:4)
    benchmark = [:unconstrianed, :constrianed, :application]
    basealgorithm =  [:NSGA2, :CCMO, :SPEA2, :SMS_EMOA]
    fcalls = 100_000
    N = 100

    allparams = @strdict seed fnum benchmark basealgorithm fcalls N

    dicts = dict_list(allparams)
    for d = dicts
        @info "Running..."
        display(d)
        algorithm, f, bounds = get_parms(d)
        res = run_algorithm(algorithm, f, bounds)
        archive = algorithm.parameters.archive

        if isempty(archive)
            @warn "Archive is empty (unable satisfy preferences)."
            continue
        end

        plt = plot_res(archive, res, d)

        @info "Saving data and Plot..."

        results = Dict("F" => fvals(archive),
                       "X" => positions(archive),
                       "fmin" => ideal(res.population),
                       "fmax" => nadir(res.population),
                       "benchmark" => d["benchmark"],
                       "fnum" => d["fnum"],
                       "basealgorithm" => d["basealgorithm"],
                       "HVROI" => compute_roi_hv(archive, d),
                       "seed" => d["seed"]
                      )
        wsave(datadir("simulations", savename(d, "bson")), results)
        wsave(plotsdir("simulations/"*string(d["basealgorithm"]), savename(d, "pdf")), plt) 
        @info "Data saved."
    end
end

main()

