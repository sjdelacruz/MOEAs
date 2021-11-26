function update_roi_archiving!(archive, population, ref_points, weight_points, δ_r, δ_w; max_archive_size=2length(population))
    if isempty(population) || (isempty(ref_points) && isempty(weight_points))
        # nothing to do
        return
    end
    
    append!(archive, population)
    mask = Metaheuristics.get_non_dominated_solutions_perm(archive)

    if isempty(mask)
        empty!(archive)
    end

    fmin = ideal(fvals(archive))
    fmax = nadir(fvals(archive))

    d = cosine_dist
    next = zeros(Bool, length(archive))
    r = ref_points
    w = weight_points
    for i in mask
        s = archive[i]

        # filtering based on ref points
        close_to_ref = any( i -> d(r[i], fval(s)) <= δ_r[i], eachindex(r))
        if close_to_ref
            next[i] = true
            continue
        end

        # filtering based on ref points
        close_to_weight = any( i -> d(w[i], (fval(s) - fmin) ./ (fmax - fmin)) <= δ_w[i], eachindex(w))

        if close_to_weight
            next[i] = true
        end
    end

    deleteat!(archive, .!next) 
    if length(archive) < max_archive_size
        return
    end
    
    # SPEA2 truncation
    del  = Metaheuristics.truncation(archive, length(archive) - max_archive_size)
    deleteat!(archive, del) 

end

