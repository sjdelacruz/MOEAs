using Metaheuristics
using Random
include("crowding-distance.jl")

mutable struct CCMO_NSGAII <: Metaheuristics.AbstractParameters
    fhelper::Metaheuristics.AbstractProblem
    M::Int
    N::Int
    η_cr::Float64
    p_cr::Float64
    η_m::Float64
    p_m::Float64
    phelper::Vector{Metaheuristics.xFgh_indiv}

end

function CCMO_NSGAII(M,fhelper_;N = 100,
    η_cr = 20,
    p_cr = 0.9,
    η_m = 20,
    p_m = -1,
    information = Information(),
    options = Options(),)

    fhelper = Problem(fhelper_[1], Array(fhelper_[2]))

    parameters = CCMO_NSGAII(fhelper,M,N, η_cr, p_cr, η_m, p_m, [])

    alg = Metaheuristics.Algorithm(
        parameters,
        information = information,
        options = options,)
    
    alg
end

function Metaheuristics.initialize!(
    status,
    parameters::CCMO_NSGAII,
    problem::Metaheuristics.AbstractProblem,
    information::Information,
    options::Options,
    args...;
    kargs...
)
    D = size(problem.bounds,2)

    if parameters.p_m < 0.0
        parameters.p_m = 1.0 / D
    end

    if options.iterations == 0
        options.iterations = 500
    end

    if options.f_calls_limit == 0
        options.f_calls_limit = options.iterations * parameters.N + 1
    end

    status = Metaheuristics.gen_initial_state(problem,parameters,information,options,status)

    

    #Generating population based on the helper problem (problem without constraints)
    parameters.phelper = Metaheuristics.generate_population(parameters.N, parameters.fhelper,ε=options.h_tol)

    status
end

function Metaheuristics.update_state!(
    status,
    parameters::CCMO_NSGAII,
    problem,
    information,
    options,
    args...;
    kargs...)

    phelper = parameters.phelper
    N = parameters.N
    middle = floor(Int, N/2)

    #For the first subset
    I = randperm(N)[1:middle]
    J = randperm(N)[1:middle]

    #For the second subset
    K = randperm(N)[1:middle]
    L = randperm(N)[1:middle]

    #Structure to save offspring
    Off1 = empty(status.population)
    Off2 = empty(status.population)

    for i = 1:2:middle

        #Select two solutions for the first subset and generate offspring
        p1a = Metaheuristics.tournament_selection(status.population, I[i])
        p1b = Metaheuristics.tournament_selection(status.population, J[i])
        p1_offspring1, p1_offspring2 = reproduction(p1a, p1b, parameters, problem)
       
        # save offsprings of forigin
        push!(Off1, p1_offspring1, p1_offspring2)

        #Select two solutions for the second subset and generate offspring
        pc = Metaheuristics.tournament_selection(phelper, K[i])
        pd = Metaheuristics.tournament_selection(phelper, L[i])
        p2_offspring1, p2_offspring2 = reproduction(pc, pd, parameters, parameters.fhelper)
       
        # save offsprings of forigin
        push!(Off2, p2_offspring1, p2_offspring2)
    end

    #Weak cooperation
    status.population = vcat(status.population, Off1,Off2)
    parameters.phelper = vcat(phelper,Off1,Off2)
    
    #Updating
    environmental_selection!(status.population, parameters)
    environmental_selection!(parameters.phelper, parameters)
end
    
function stop_criteria_ctaea(
    status,
    parameters::CCMO_NSGAII,
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
    parameters::CCMO_NSGAII,
    problem,
    information,
    options,
    args...;
    kargs...
)
    
    status.final_time = time()
end

function reproduction(pa, pb, parameters::CCMO_NSGAII, problem)
    
    c1, c2 = Metaheuristics.GA_reproduction(get_position(pa),
                             get_position(pb),
                             problem.bounds;
                             η_cr = parameters.η_cr,
                             p_cr = parameters.p_cr,
                             η_m = parameters.η_m,
                             p_m = parameters.p_m)

    Metaheuristics.create_solution(c1, problem), Metaheuristics.create_solution(c2, problem) 
end

function environmental_selection!(population, parameters::CCMO_NSGAII)
    truncate_population!(population, parameters.N)
end


