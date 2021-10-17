
using Distances
using Metaheuristics

include("HardTestProblems/HardTestProblems.jl")
using .HardTestProblems

include("CTAEA.jl")

function main()
    
    #Objective function and a set of parameters
    f, conf = get_RW_MOP_problem(1)
    fx, gx, hx = f(conf[:xmin])
    vio = sum(max.(0.0, gx .<= 0)) + sum(max.(0.0, abs.(hx) .- 1e-2 .<= 0))
    
    #Dimensions
    D = conf[:n]

    #Setting bounds
    xmin, xmax = conf[:xmin],conf[:xmax]
    _bounds = (xmin, xmax)
    bounds = vcat(transpose.(_bounds)...)

    #Information - doubt
    information = Information(f_optimum = 0.0)

    options = Options(f_calls_limit = 9000*10, f_tol = 1e-5)

    # reference points (Das and Dennis's method)
    weights = gen_ref_dirs(D, 10)

    # algoritmo a utilizar
    algorithm = CTAEA(weights,information = information, options = options)

    # Start optimization process
    resultado = optimize(f, bounds, algorithm)

end

main()