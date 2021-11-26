mutable struct ROIs{T} <: Metaheuristics.AbstractParameters
    parameters::T
    weight_points::Vector{Vector{Float64}}
    ref_points::Vector{Vector{Float64}}
    δ_r::Vector{Float64}
    δ_w::Vector{Float64}
    archive::Vector{Metaheuristics.xFgh_indiv}
    max_archive_size::Int
end

function ROIs(algorithm = CCMO_NSGAII();
        weight_points=Vector{Float64}[],
        ref_points=Vector{Float64}[],
        δ_r = zeros(0),
        δ_w = zeros(0),
        max_archive_size = 0
    )

    length(ref_points) != length(δ_r) && error("|ref_points| is different to |δ_r|")
    length(weight_points) != length(δ_w) && error("|weight_points| is different to |δ_w|")

    parameters = ROIs(algorithm.parameters, weight_points, ref_points, δ_r, δ_w, Metaheuristics.xFgh_indiv[],max_archive_size)


    return Metaheuristics.Algorithm(parameters, information = algorithm.information, options = algorithm.options)

end


function Metaheuristics.initialize!(
        status,
        parameters::ROIs,
        args...;
        kargs...
    )

    st = Metaheuristics.initialize!(status, parameters.parameters,args...;kargs...)
    if parameters.max_archive_size == 0
        parameters.max_archive_size = 2parameters.parameters.N
    end

    st
    
end

function Metaheuristics.update_state!(
        status,
        parameters::ROIs,
        args...;
        kargs...
    )

    Metaheuristics.update_state!(status, parameters.parameters,args...;kargs...)

    archive = parameters.archive
    population = status.population

    ref_points = parameters.ref_points
    δ_r = parameters.δ_r

    weight_points = parameters.weight_points
    δ_w = parameters.δ_w
    max_archive_size = parameters.max_archive_size

    update_roi_archiving!(archive, population, ref_points, weight_points, δ_r, δ_w;max_archive_size = max_archive_size)
end

function Metaheuristics.final_stage!(
        status,
        parameters::ROIs,
        args...;
        kargs...
    )

    Metaheuristics.final_stage!(status, parameters.parameters,args...;kargs...)
end
