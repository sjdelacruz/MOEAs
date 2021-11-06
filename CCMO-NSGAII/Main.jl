#Name = Sebastián José de la Cruz Martínez
#Main function to call CCMO_NSGA2 using metaheuristcs package
#Problem with C1-DTLZ3 AND DTLZ3, dimensions D in the paper are wrong,
# not match with the parameters defined in the suite.

try
    using HardTestProblems
    using Distances
    using Plots
    pyplot()
    using HardTestProblems
    using Metaheuristics
    

catch
    import Pkg; Pkg.add("HardTestProblems")    
    Pkg.rm("Metaheuristics")
    Pkg.add(url="https://github.com/jmejia8/Metaheuristics.jl.git#develop")
    Pkg.add("Plots")
    Pkg.add("PyPlot")
end

include("CCMO_NSGAII.jl")

function main() 

    M=3;
    f, bounds, front = Metaheuristics.TestProblems.C1_DTLZ3(M)
    fhelper = Metaheuristics.TestProblems.DTLZ3(M)
    display(front)
    #plt3d= Plots.plot(fs[:,1],fs[:,2], fs[:,3], seriestype=:scatter, markersize = 7,title = "Pareto Front")
    #savefig("Pf.png")

    options = Options(f_calls_limit=100000, debug=false)
    algorithm = CCMO_NSGAII(M,fhelper,N = 105, p_cr = 0.85, options = options)

    # Start optimization process
    resultado = optimize(f, bounds, algorithm)
    display(resultado)
end

main()
