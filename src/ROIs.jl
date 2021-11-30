mutable struct ROIs{T} <: Metaheuristics.AbstractParameters
    parameters::T
    weight_points::Vector{Vector{Float64}}
    δ_w::Vector{Float64}
    archive::Vector{Metaheuristics.xFgh_indiv}
    max_archive_size::Int
end

function ROIs(algorithm = CCMO_NSGAII();
        weight_points=Vector{Float64}[],
        δ_w = zeros(0),
        max_archive_size = 0
    )

    length(weight_points) != length(δ_w) && error("|weight_points| is different to |δ_w|")


    δ_w /= 2
    δ_w = 1 .- cos.(π/2*δ_w)


    parameters = ROIs(algorithm.parameters, weight_points, δ_w, Metaheuristics.xFgh_indiv[],max_archive_size)


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
        parameters.max_archive_size = parameters.parameters.N
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


    weight_points = parameters.weight_points
    δ_w = parameters.δ_w
    max_archive_size = parameters.max_archive_size

    update_roi_archiving!(archive, population, weight_points, δ_w;max_archive_size = max_archive_size)
end

function Metaheuristics.final_stage!(
        status,
        parameters::ROIs,
        args...;
        kargs...
    )

    Metaheuristics.final_stage!(status, parameters.parameters,args...;kargs...)

    # updating archive, remove those solutions not in ROI
    archive = parameters.archive
    if isempty(archive)
        return
    end

    w = parameters.weight_points
    δ_w = parameters.δ_w

    fmin = ideal(status.population)
    fmax = nadir(status.population)

    g_roi = compute_rio_vio(archive, fmin, fmax, w, δ_w)
    mask = g_roi .!= 0
    deleteat!(archive, mask)
end
