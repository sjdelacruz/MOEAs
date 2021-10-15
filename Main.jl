#import Pkg; 
#Pkg.add("Distances")
using Distances
using Metaheuristics
include("Multiobjective/RW_MOP_2021/RW_MOP_2021.jl")
include("CTAEA.jl")


function main()

    # objective function
    f, bounds, pareto_solutions = TestProblems.get_problem(:ZDT3);
    nobjectives = 2
    npartitions = 100
    nrestricciones = 0

    # reference points (Das and Dennis's method)
    Wv = gen_ref_dirs(nobjectives, npartitions)

    # algoritmo a utilizar
    ctaea = CTAEA(Wv, options = Options(debug=false, iterations= 250))

    # inicializar el proceso de optimización
    resultado = optimize(f, bounds, ctaea)

    #Mostrar resultados
    display(status_moead)
end

main()