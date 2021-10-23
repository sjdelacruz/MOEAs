using Metaheuristics
include("UpdateDA.jl")
include("UpdateCA.jl")

mutable struct CTAEA <: Metaheuristics.AbstractParameters
    nobjectives::Int
    N::Int #individuos
    weights
end

function CTAEA(weights; N = 100, information = Information(), options = Options())
    
    if isempty(weights)
        error("Provide weighted vectors")
    end

    nobjectives = length(weights[1])

    parameters = CTAEA(nobjectives, N, weights)

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

    return Metaheuristics.gen_initial_state(problem,parameters,information,options,status)

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

    #Updating archives
    #At the beginning CA is empty
    CA = UpdateCA(empty(population), population, weights, N);
    # DA = UpdateDA(CA, empty(population), population, weights); 


    # remove following lines when this function works correctly
    @info "I'm in update_state :)"
    status.stop = true
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
