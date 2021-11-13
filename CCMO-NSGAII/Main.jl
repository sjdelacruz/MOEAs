#Name = Sebastián José de la Cruz Martínez
#Main function to call CCMO_NSGA2 using metaheuristcs package
#Problem with C1-DTLZ3 AND DTLZ3, dimensions D in the paper are wrong,
# not match with the parameters defined in the suite.

try
    using HardTestProblems
catch
    import Pkg; Pkg.add("HardTestProblems")    
    Pkg.rm("Metaheuristics")
    Pkg.add(url="https://github.com/jmejia8/Metaheuristics.jl.git#develop")
    Pkg.add("Plots")
    Pkg.add("PyPlot")
end

using Distances
using Plots
pyplot()
using HardTestProblems
using Metaheuristics
using Statistics
using LinearAlgebra

include("CCMO_NSGAII.jl")

function test_preferences()
    M = 3
    f, bounds, front = Metaheuristics.TestProblems.DTLZ1(M);
    D = size(bounds,2)

    # preferences defined as weights
    preferences = [[0.2,0.7, 0.1],
                   [0.25, 0.25, 0.5],
                   [0.7, 0.2, 0.1],
                  ]
    # the threshold (cosine distance)
    δ = fill(0.01, length(preferences))

    #= preferences as constraints
    f(x,w=preferences,δ=δ) = begin
        fx = ff(x)[1]
        g = [cosine_dist(fx, ww) - δ for ww in w] 
        return fx, [minimum(g)], [0.0]
    end
    =#

    options = Options(f_calls_limit=100000,iterations=10000, debug=false)
    algorithm = CCMO_NSGAII(;N = 105, p_cr = 0.9, p_m= (1.0/D), options = options, preferences,δ)
    res = Metaheuristics.optimize(f, bounds, algorithm)
    display(res)
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
end


function main() 

    M=3;
    D=7
    f, bounds, front = Metaheuristics.TestProblems.C1_DTLZ3(M)
    bounds = bounds[:,1:D]
    
    #plt3d= Plots.plot(fs[:,1],fs[:,2], fs[:,3], seriestype=:scatter, markersize = 7,title = "Pareto Front")
    #savefig("Pf.png")

    igd_values = []

    for i in 1:30
        options = Options(f_calls_limit=100000, debug=false)
        algorithm = CCMO_NSGAII(N = 105, p_cr = 1.0, p_m= (1.0/D), options = options)
        resultado = Metaheuristics.optimize(f, bounds, algorithm)
        igd_local = Metaheuristics.PerformanceIndicators.igd(resultado,front)
        push!(igd_values, igd_local)
    end

    igd_mean = mean(igd_values)
    println("#Mean IGD: ", igd_mean)
end

# main()
test_preferences()
