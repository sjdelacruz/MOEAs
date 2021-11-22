using DrWatson
using Metaheuristics
@quickactivate "MOEAs"


include("../src/MOEAs.jl")

function run_algorithm(algorithm, f, bounds)
    res = Metaheuristics.optimize(f, bounds, algorithm)
    display(res)
    
end


function main()
    M = 3
    f, bounds, front = Metaheuristics.TestProblems.C2_DTLZ2(M);
    D = size(bounds,2)

    # preferences defined as weights
    preferences = [[0.2,0.7, 0.1],
                   [0.25, 0.25, 0.5],
                   [0.7, 0.2, 0.1],
                  ]
    # the threshold (cosine distance)
    δ = fill(0.01, length(preferences))


    options = Options(f_calls_limit=100000,iterations=10000, debug=false)
    algorithm = MOEAs.CCMO_NSGAII(;N = 105, p_cr = 0.9, p_m= (1.0/D), options = options, preferences,δ)
    run_algorithm(algorithm, f, bounds)
    #=
    p = plot(layout=(1,2),xlabel="f₁", ylabel="f₂", zlabel="f₃",xlim=[0,1],ylim=[0,1], zlim=[0,1])
    fs = fvals(front)
    wireframe!(p[1], fs[:,1], fs[:,2], fs[:,3], linecolor=:lightgray, fillalpha=0)
    # wireframe!(p[2], fs[:,1], fs[:,2], fs[:,3], linecolor=:lightgray, fillalpha=0)
    fs = fvals(res)
    scatter!(p[1], fs[:,1], fs[:,2], fs[:,3], label="Result")
    scatter!(p[2], fs[:,1], fs[:,2], fs[:,3], label="Result")

    for (i,w) in enumerate(preferences)
        t = range(0,1, length=50)
        ww = w / norm(w)
        line = zeros(3)' .+ t.*ww'
        plot!(p[1], line[:,1], line[:,2], line[:,3], label="Preference $i", color=:blue, lw=2)
        plot!(p[2], line[:,1], line[:,2], line[:,3], label="Preference $i", color=:blue, lw=2)
    end
    p
    =#
end

main()

