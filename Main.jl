
using Distances
using Metaheuristics

try
    using HardTestProblems
catch
    # Install module
    import Pkg; Pkg.add(url="https://github.com/jmejia8/HardTestProblems.jl")
end
using HardTestProblems

include("CTAEA.jl")

function main()
    
    #Objective function and a set of parameters
    f, conf = get_RW_MOP_problem(1)
    fx, gx, hx = f(conf[:xmin])
    vio = sum(max.(0.0, gx .<= 0)) + sum(max.(0.0, abs.(hx) .- 1e-2 .<= 0))
    
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
    algorithm = CTAEA(weights, options = options)

    # Start optimization process
    resultado = optimize(f, bounds, algorithm)

end

main()
