using Distances
using Metaheuristics
using Plots;pyplot()

try
    using HardTestProblems
catch
    # Install module
    import Pkg; Pkg.add("HardTestProblems")
end
using HardTestProblems

include("CTAEA.jl")

function main_old()
    
    #Objective function and a set of parameters
    f, conf = get_RW_MOP_problem(3)

    
    #Dimensions
    D = conf[:n]
    
    # objectives
    M = conf[:fn]

    #Setting bounds
    xmin, xmax = conf[:xmin],conf[:xmax]
    bounds = [xmin xmax]

    options = Options(f_calls_limit = 9000*10, debug=false)

    # reference points (Das and Dennis's method)
    weights = gen_ref_dirs(M, 100)

    # algoritmo a utilizar
    algorithm = CTAEA(N = 100, weights, options = options)

    # Start optimization process
    resultado = optimize(f, bounds, algorithm)
    #fvals(resultado)

end

function main() 

    f, bounds, front = Metaheuristics.TestProblems.C1_DTLZ3()

    M = 3
    options = Options(iterations=1000, debug=false)

    # reference points (Das and Dennis's method)
    weights = gen_ref_dirs(M, 12)

    # algoritmo a utilizar
    algorithm = CTAEA(weights, N = length(weights), η_cr = 30, η_m = 20, options = options)

    # Start optimization process
    resultado = optimize(f, bounds, algorithm)
    display(resultado)


    fs = fvals(front)
    wireframe(fs[:,1], fs[:,2], fs[:,3])
    fs = fvals(resultado)
    scatter!(fs[:,1], fs[:,2], fs[:,3])
    #fvals(resultado)
end


main()
