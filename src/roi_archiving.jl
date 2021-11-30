using LinearAlgebra

function update_roi_archiving!(archive, population, weight_points, δ_w; max_archive_size=2length(population))
    if isempty(population) ||  isempty(weight_points)
        # nothing to do
        return
    end

    feasible_sols = Metaheuristics.is_feasible.(population)
    if !any(feasible_sols)
        # we cannot work on infeasible solutions
        return
    end

    non_dominated = Metaheuristics.get_non_dominated_solutions_perm(population)
    if length(non_dominated) < length(population)
        # the front is not well distributed
        return
    end
    
    #empty!(archive)

    fmin = ideal(population)
    fmax = nadir(population)


    w = weight_points


    fs = fvals(population)
    M = size(fs,2)
    extremas = vcat([argmin(fs[:,i]) for i in 1:M ], [argmax(fs[:,i]) for i in 1:M ])
    unique!(extremas)

    g_roi = compute_rio_vio(population, fmin, fmax, weight_points, δ_w)
    for (l,s) in enumerate(population)
        # filtering based on weight_points
        gx = g_roi[l]
        gx <= 0 && push!(archive, s)

        # no manipulate extrema points
        if l in extremas
            continue
        end

        s.sum_violations += max(gx, 0)
        s.is_feasible = s.sum_violations <= 0


    end

    truncate_archive!(archive,  weight_points, δ_w; max_archive_size)
end

function compute_rio_vio(population, fmin, fmax, w, δ_w)

    d = cosine_dist
    g = zeros(length(population))
    for (l,s) in enumerate(population)


        # filtering based on weight_points
        gx = minimum( i -> d(w[i], (fval(s) - fmin) ./ (fmax - fmin)) - δ_w[i], eachindex(w))
        g[l] = max(gx,0)

    end
    return g
end


function truncate_archive!(archive,  weight_points, δ_w; max_archive_size=length(population))
    unique!(archive)
    mask = Metaheuristics.get_non_dominated_solutions_perm(archive)

    if isempty(mask)
        empty!(archive)
        return
    elseif length(mask) <= max_archive_size
        return
    end

    next = zeros(Bool, length(archive))
    next[mask] .= true
    deleteat!(archive, .!next)
    
    if length(archive) > max_archive_size
        # SPEA2 truncation/ knn
        del  = Metaheuristics.truncation(archive, length(archive) - max_archive_size)
        deleteat!(archive, del) 
        return
    end
end

