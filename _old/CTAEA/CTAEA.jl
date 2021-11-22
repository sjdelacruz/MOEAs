using Metaheuristics
include("UpdateDA.jl")
include("UpdateCA.jl")
include("MatingSelection.jl")

mutable struct CTAEA <: Metaheuristics.AbstractParameters
    nobjectives::Int
    N::Int #individuos
    # crossover
    η_cr::Float64
    p_cr::Float64
    η_m::Float64
    p_m::Float64
    weights
    CA::Vector{Metaheuristics.xFgh_indiv}
    DA::Vector{Metaheuristics.xFgh_indiv}
end

function CTAEA(weights; N = 100, η_cr = 20, p_cr=0.9, η_m=15, p_m=0.1, information = Information(), options = Options())
    
    if isempty(weights)
        error("Provide weighted vectors")
    end

    nobjectives = length(weights[1])

    parameters = CTAEA(nobjectives, N, η_cr, p_cr, η_m, p_m, weights, [], [])

    alg = Metaheuristics.Algorithm(
        parameters,
        information = information,
        options = options,)
    
    alg
end

function Metaheuristics.initialize!(
    status,
    parameters::CTAEA,
    problem::Metaheuristics.AbstractProblem,
    information::Information,
    options::Options,
    args...;
    kargs...
)
        
    if options.iterations == 0
        options.iterations = 500
    end

    if options.f_calls_limit == 0
        options.f_calls_limit = options.iterations * parameters.N + 1
    end


    options.debug && @info "Initializing Archives."


    status = Metaheuristics.gen_initial_state(problem,parameters,information,options,status)

    population = status.population
    weights = parameters.weights
    N = parameters.N
    #Updating archives
    #At the beginning CA is empty
    parameters.CA = UpdateCA(empty(population), population, weights, N)
    parameters.DA = UpdateDA(parameters.CA, empty(population), population, weights)

    status
end

function Metaheuristics.update_state!(
    status,
    parameters::CTAEA,
    problem,
    information,
    options,
    args...;
    kargs...)

    N = parameters.N
    D = size(problem.bounds, 2)
    population = status.population
    weights = parameters.weights
    CA = parameters.CA
    DA = parameters.DA

    ## mating pool choosing
    # calculate the ratio of non-dominated solutions of CA and DA in Hm
    Hm = vcat(CA,DA)                         
    # [FrontNo,~]=NDSort(Hm.objs,inf);
    Metaheuristics.fast_non_dominated_sort!(Hm)
    FrontNo = map(s -> s.rank, Hm)
    FrontNo_C=FrontNo[1:ceil(Int, length(Hm)/2)]

    Nc = length(findall(FrontNo_C .==1 ))
    Pc = Nc/length(Hm)
    FrontNo_D = FrontNo[ceil(Int, length(Hm)/2)+1:length(Hm)]
    Nd = length(findall(FrontNo_D .==1 ))
    Pd = Nd/length(Hm)

    # calculate the proportion of non-dominated solutions in CA
    Metaheuristics.fast_non_dominated_sort!(CA)
    NC = length(findall(FrontNo .== 1))
    # PC denotes the proportion of non-dominated solutions in CA,it is different from Pc
    PC = NC/length(CA)

    #reproduction
    Q = empty(population)


    W = parameters.weights
    for i in 1:length(W)
        if Pc > Pd
            P1 = MatingSelection(CA)
        else
            P1 = MatingSelection(DA)
        end

        pf = rand();
        if pf < PC
            P2 = MatingSelection(CA);
        else
            P2 = MatingSelection(DA);
        end

        x1 = Metaheuristics.get_position(P1)
        x2 = Metaheuristics.get_position(P2)

        c, _ = Metaheuristics.SBX_crossover(x1, x2, problem.bounds, parameters.η_cr, parameters.p_cr)
        Metaheuristics.polynomial_mutation!(c,problem.bounds,parameters.η_m, parameters.p_m)
        Metaheuristics.reset_to_violated_bounds!(c, problem.bounds)
        offspring = Metaheuristics.create_solution(c, problem)

        push!(Q, offspring)
    end

    parameters.CA = UpdateCA(parameters.CA, Q, weights, N)
    parameters.DA = UpdateDA(parameters.CA, parameters.DA, Q, weights)


    # remove solutions with NAN or Inf values
    status.population = parameters.CA
    mask = map(s -> any(.!isfinite.(fval(s)) .| isnan.(fval(s))), status.population)
    status.population = status.population[.!mask]

    #status.stop = true
end
    


function stop_criteria_ctaea(
    status,
    parameters::CTAEA,
    problem,
    information,
    options,
    args...;
    kargs...
    )

    return status.iteration > options.iterations
end

function Metaheuristics.final_stage!(
    status,
    parameters::CTAEA,
    problem,
    information,
    options,
    args...;
    kargs...
)
    
    status.final_time = time()
end
