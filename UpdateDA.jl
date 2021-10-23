using Metaheuristics
using Distances
# DA=UpdateDA(CA,[],status.population,parameters.Wv);
#
"""
    UpdateDA(CA,DA,Q,W)

Update archive where CA, DA are two populations, W contains the weights.
The output is like a population population.

**Warning:** This implementation is like the matlab code, i.e., slow. Performance update is required.
**Warning:** rank are updated here, be aware.
"""
function UpdateDA(CA,DA,Q,W_)

    W = hcat(W_...)


    # S is the set used for output
    S=empty(Q)
    Hd= vcat(DA,Q)
    N=size(W,1)


    Region_Hd_ = argmax(1 .- pairwise(cosine_dist, fvals(Hd)', W), dims = 2)
    Region_Hd = map(c -> c.I[2], Region_Hd_[:,1])
    # associat the individuals in Hd and CA with their corresponding subregions
    Region_CA_ = argmax(1 .- pairwise(cosine_dist, fvals(CA)', W), dims = 2)
    Region_CA = map(c -> c.I[2], Region_CA_[:,1])


    itr = 1

    while length(S)<N

        # here i denotes the order of the subregion
        for i in 1:N              
            
            # current_c denotes that the current_c_th individual is/are already in the ith region
            current_c=findall(Region_CA .== i);                       
            if length(current_c)<itr 
                # j denotes the number of solutions from Hd that need to join into the region(i)
                for j in 1:itr-length(current_c)                   
                    current_d = findall(Region_Hd .== i);
                    if isempty(current_d)
                        break
                    end

                    #[FrontNO,~]=NDSort(Hd(current_d).objs,inf);
                    Metaheuristics.fast_non_dominated_sort!(Hd[current_d])
                    FrontNO = map(s -> s.rank, Hd[current_d])

                    #O is the set of nondominated solutions from region(i) in Hd
                    O=Hd[current_d[FrontNO .==1]];            
                    Z = ideal(O)

                    Region_O_ = argmax(1 .-pairwise(cosine_dist, fvals(O)',W), dims = 2)
                    Region_O = map(c -> c.I[2], Region_O_[:,1])

                    g_tch=maximum(abs.(fvals(O) .- Z') ./ W[:, Region_O]', dims = 2)[:,1]
                    order=argmin(g_tch)
                    x_best=O[order]


                    mask = findall(s -> s === x_best, Hd[current_d])
                    deleteat!(Hd, mask)

                    if !isempty(Hd)
                        Region_Hd_ = argmax(1 .-pairwise(cosine_dist, fvals(Hd)',W), dims =2)
                        Region_Hd = map(c -> c.I[2], Region_Hd_[:,1])
                    end    

                    # add the best individual into S
                    if length(S)<N                          
                        push!(S, x_best)
                    end
                end
            end

            if length(S)==N
                break
            end
        end
        itr += 1;
    end

    return S
end

#=
function test()
    f, bounds, pareto_solutions = Metaheuristics.TestProblems.get_problem(:ZDT3);
    nobjectives = 2
    npartitions = 100

    Wv = gen_ref_dirs(nobjectives, npartitions)
    @time DA=UpdateDA(rand(pareto_solutions, 10), pareto_solutions, empty(pareto_solutions) ,Wv);

end

test()
=#
