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

include("CCMO_NSGAII.jl")

function main() 

    M=3;
    D=7
    f, bounds, front = Metaheuristics.TestProblems.C1_DTLZ3(M)
    bounds = bounds[:,1:D]
    
    fhelper = Metaheuristics.TestProblems.DTLZ3(M)
    #plt3d= Plots.plot(fs[:,1],fs[:,2], fs[:,3], seriestype=:scatter, markersize = 7,title = "Pareto Front")
    #savefig("Pf.png")

    igd_values = []

    for i in 1:30
        options = Options(f_calls_limit=100000, debug=false)
        algorithm = CCMO_NSGAII(M,D,fhelper, N = 105, p_cr = 1.0, p_m= (1.0/D), options = options)
        resultado = Metaheuristics.optimize(f, bounds, algorithm)
        igd_local = Metaheuristics.PerformanceIndicators.igd(resultado,front)
        push!(igd_values, igd_local)
    end

    igd_mean = mean(igd_values)
    println("#Mean IGD: ", igd_mean)
end

main()
