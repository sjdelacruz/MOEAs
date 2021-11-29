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


    d = cosine_dist
    dd = norm(fmin - fmax)
    next = zeros(Bool, length(archive))
    w = weight_points


    fs = fvals(population)
    M = size(fs,2)
    extremas = vcat([argmin(fs[:,i]) for i in 1:M ], [argmax(fs[:,i]) for i in 1:M ])
    unique!(extremas)
    for (l,s) in enumerate(population)
        if l in extremas
            continue
        end


        # filtering based on weight_points
        gx = minimum( i -> d(w[i], (fval(s) - fmin) ./ (fmax - fmin)) - δ_w[i], eachindex(w))
        gx <= 0 && push!(archive, s)
        s.sum_violations += max(gx, 0)
        s.is_feasible = s.sum_violations <= 0


    end

    truncate_archive(archive,  weight_points, δ_w; max_archive_size)
end

function truncate_archive(archive,  weight_points, δ_w; max_archive_size=2length(population))
    unique!(archive)
    mask = Metaheuristics.get_non_dominated_solutions_perm(archive)

    if isempty(mask)
        empty!(archive)
        return
    elseif length(mask) == max_archive_size
        return
    end

    next = zeros(Bool, length(archive))
    next[mask] .= true
    deleteat!(archive, .!next)
    
    if length(archive) > max_archive_size
        # SPEA2 truncation
        del  = Metaheuristics.truncation(archive, length(archive) - max_archive_size)
        deleteat!(archive, del) 
        return
    end
end

