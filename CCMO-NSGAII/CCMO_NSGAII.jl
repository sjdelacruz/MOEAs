using Metaheuristics
using Random

mutable struct CCMO_NSGAII <: Metaheuristics.AbstractNSGA
    N::Int
    η_cr::Float64
    p_cr::Float64
    η_m::Float64
    p_m::Float64
    # helper population
    phelper::Vector{Metaheuristics.xFgh_indiv}

end

function CCMO_NSGAII(;
        N = 100,
        η_cr = 20,
        p_cr = 0.9,
        η_m = 20,
        p_m = -1,
        information = Information(),
        options = Options(),
    )

    parameters = CCMO_NSGAII(N, η_cr, p_cr, η_m, p_m, [])

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
    parameters.phelper = Metaheuristics.generate_population(parameters.N, problem,ε=options.h_tol)

    status
end

function Metaheuristics.update_state!(
    status::State,
    parameters::CCMO_NSGAII,
    problem::Metaheuristics.AbstractProblem,
    information::Information,
    options::Options,
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
        p1_offspring1, p1_offspring2 = Metaheuristics.reproduction(p1a, p1b, parameters, problem)
       
        # save offsprings of forigin
        push!(Off1, p1_offspring1, p1_offspring2)

        #Select two solutions for the second subset and generate offspring
        pc = Metaheuristics.tournament_selection(phelper, K[i])
        pd = Metaheuristics.tournament_selection(phelper, L[i])
        p2_offspring1, p2_offspring2 = Metaheuristics.reproduction(pc, pd, parameters, problem)
       
        # save offsprings of forigin
        push!(Off2, p2_offspring1, p2_offspring2)
    end

    # Weak cooperation
    status.population = vcat(status.population, Off1, Off2)
    # copy Off1 and Off2 (copy objects not pointers)
    parameters.phelper = vcat(phelper, deepcopy(Off1), deepcopy(Off2))
    
    #Updating
    environmental_selection!(status.population, parameters)
    environmental_selection!(parameters.phelper, parameters, false)
end


function Metaheuristics.final_stage!(
    status::State,
    parameters::CCMO_NSGAII,
    problem::Metaheuristics.AbstractProblem,
    information::Information,
    options::Options,
    args...;
    kargs...
)
    
    status.final_time = time()
end


function environmental_selection!(population, parameters::CCMO_NSGAII, consider_constrints=true)
    if consider_constrints
        # using constrained non-dominated sorting here
        Metaheuristics.truncate_population!(population, parameters.N)
        return
    end

    # using non-dominated (without constraints)
    CV = [s.sum_violations for s in parameters.phelper]
    # putting CV = 0 to ingnore constraints in helper population
    for s in parameters.phelper
        s.sum_violations = 0
    end

    Metaheuristics.truncate_population!(parameters.phelper, parameters.N)

    # restore constraint violation values
    for (i,s) in enumerate(parameters.phelper)
        s.sum_violations = CV[i]
    end
    
end


