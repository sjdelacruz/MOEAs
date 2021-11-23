using Metaheuristics
using Distances
using Random

mutable struct CCMO_NSGAII <: Metaheuristics.AbstractNSGA
    N::Int
    η_cr::Float64
    p_cr::Float64
    η_m::Float64
    p_m::Float64
    # helper population
    phelper::Vector{Metaheuristics.xFgh_indiv}
    preferences::Vector{Vector{Real}}
    δ::Vector{Float64}
    fitness::Vector{Float64}
    fitness_helper::Vector{Float64}

end

function CCMO_NSGAII(;
        N = 100,
        η_cr = 20,
        p_cr = 0.9,
        η_m = 20,
        p_m = -1,
        preferences=[],
        δ=[],
        information = Information(),
        options = Options(),
    )

    parameters = CCMO_NSGAII(N, η_cr, p_cr, η_m, p_m, [],preferences,δ,[],[])

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

    # used to compute fitness
    environmental_selection!(status.population, parameters)
    environmental_selection!(parameters.phelper, parameters, false)

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
    fitness1 = parameters.fitness
    fitness2 = parameters.fitness_helper
    N = parameters.N
    middle = floor(Int, N/2)


    #Structure to save offspring
    Off1 = empty(status.population)
    Off2 = empty(status.population)

    for i = 1:2:middle

        #Select two solutions for the first subset and generate offspring
        p1a = Metaheuristics.binary_tournament(status.population, fitness1)
        p1b = Metaheuristics.binary_tournament(status.population, fitness1)
        p1_offspring1, p1_offspring2 = Metaheuristics.reproduction(p1a, p1b, parameters, problem)
       
        # save offsprings of forigin
        push!(Off1, p1_offspring1, p1_offspring2)

        #Select two solutions for the second subset and generate offspring
        pc = Metaheuristics.binary_tournament(phelper, fitness2)
        pd = Metaheuristics.binary_tournament(phelper, fitness2)
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
    d = cosine_dist
    ws = parameters.preferences
    δ = parameters.δ
    spea2 = SPEA2().parameters
    spea2.N = parameters.N
    if consider_constrints
        tmp = copy(population)
        CV = [s.sum_violations for s in population]
        # handling preferences
        F = fvals(population)
        fmin = ideal(F)'
        fmax = nadir(F)'
        Fnorm = F#(F .- fmin) ./ (fmax - fmin)
        for (i,s) in enumerate(population)
            gg = minimum([d(Fnorm[i,:], w)-δ[j] for (j,w) in enumerate(ws)])
            s.sum_violations += max(gg, 0)
            s.is_feasible = s.sum_violations == 0
        end

        # using constrained non-dominated sorting here
        # Metaheuristics.truncate_population!(population, parameters.N)

        Metaheuristics.environmental_selection!(population,spea2)
        parameters.fitness = spea2.fitness
        #=
        for (i,s) in enumerate(tmp)
            s.sum_violations = CV[i]
            s.is_feasible = CV[i] == 0
        end
        =#

    else

        CV = [s.sum_violations for s in parameters.phelper]
        tmp = copy(parameters.phelper)
        # using non-dominated (without constraints)
        # putting CV = 0 to ingnore constraints in helper population
        for s in parameters.phelper
            s.sum_violations = 0
        end

        Metaheuristics.environmental_selection!(parameters.phelper,spea2)
        #Metaheuristics.truncate_population!(parameters.phelper, parameters.N)
        # restore constraint violation values
        for (i,s) in enumerate(tmp)
            s.sum_violations = CV[i]
        end
        parameters.fitness_helper = spea2.fitness
    end


    
end


