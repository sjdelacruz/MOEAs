using LinearAlgebra

function update_roi_archiving!(archive, population, ref_points, weight_points, δ_r, δ_w; max_archive_size=2length(population))
    if isempty(population) || (isempty(ref_points) && isempty(weight_points))
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
    
    empty!(archive)

    fmin = ideal(population)
    fmax = nadir(population)


    d = cosine_dist
    dd = norm(fmin - fmax)
    next = zeros(Bool, length(archive))
    r = ref_points
    w = weight_points


    fs = fvals(population)
    M = size(fs,2)
    extremas = vcat([argmin(fs[:,i]) for i in 1:M ], [argmax(fs[:,i]) for i in 1:M ])
    unique!(extremas)
    for (l,s) in enumerate(population)
        if l in extremas
            continue
        end

        #= filtering based on ref points
        close_to_ref = any( i -> norm(r[i] - fval(s)) <= dd*δ_r[i], eachindex(r))
        if close_to_ref
            next[i] = true
            continue
        end
        =#

        # filtering based on ref points
        gx = minimum( i -> d(w[i], (fval(s) - fmin) ./ (fmax - fmin)) - δ_w[i], eachindex(w))
        gx < 0 && push!(archive, deepcopy(s))
        s.sum_violations += max(gx, 0)
        s.is_feasible = s.sum_violations <= 0


        #=
        if close_to_weight
            next[i] = true
        end
        =#
    end
    

    #=
    append!(archive, population)
    unique!(archive)
    #=
    mask = Metaheuristics.get_non_dominated_solutions_perm(archive)

    @show length(archive)
    if isempty(mask)
        empty!(archive)
    end
    =#


    d = cosine_dist
    dd = norm(fmin- fmax)
    next = zeros(Bool, length(archive))
    r = ref_points
    w = weight_points
    for i in mask
        s = archive[i]

        # filtering based on ref points
        close_to_ref = any( i -> norm(r[i] - fval(s)) <= dd*δ_r[i], eachindex(r))
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

    if length(archive) < 0.1max_archive_size
        return
    end
    
    @info "population"
    if length(archive) > max_archive_size
        # SPEA2 truncation
        del  = Metaheuristics.truncation(archive, length(archive) - max_archive_size)
        deleteat!(archive, del) 
        return
    end

    @info "clonning"
    @show length(archive)

    N = length(population)
    fs = fvals(population)
    M = size(fs,2)
    extr = Int[argmin(fs[:,i]) for i in 1:M ]
    extr2 = Int[argmax(fs[:,i]) for i in 1:M ]
    append!(extr, extr2)
    unique!(extr)
    ss = population[extr]
    empty!(population)
    append!(population, ss)
    
    i = 1
    while length(population) < N
        push!(population, archive[i])
        i += 1
        i = i > length(archive) ? 1 : i
    end
    =#


end

