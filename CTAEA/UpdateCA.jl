using Distances
using Metaheuristics
# using StatsBase

#N Number of solutions (at the beginning is the initial population)

function UpdateCA(CA,Q,W_,N)
    
    #Change to dimensions weighted vectors
    W = hcat(W_...)


    #Output
    S = empty(Q)

    #Hybrid population
    Hc = vcat(CA,Q)
     
    # Feasible solutions 
    Sc = Hc[findall(s -> s.is_feasible, Hc)]

    #Number of feasible solutions   
    feasibles = length(Sc)
    
    #1st condition
    if feasibles == N
        CA = Sc

        return CA
    # 2nd condition
    elseif feasibles >= N
        
        #Non dominated sorting to Sc
        Metaheuristics.fast_non_dominated_sort!(Sc)
        FrontNO = map(s -> s.rank, Sc)

        nfront = length(FrontNO)

        #Fill S with fronts i=1 ... nfront
        for i in 1:nfront 
            fp = filter(s -> s.rank == FrontNO[i], Sc)
            S = append!(S, fp)
            if length(S) >= N
                break
            end
        end

        while length(S) > N

            SPopObj = fvals(S)
            # Ideal and nadir points
            zideal = ideal(SPopObj)
            znad = nadir(SPopObj)
            
            # Normalize_points
            SPopObj = ((SPopObj .- zideal') ./ (znad - zideal)')

            # association of S and W 
            Region_S_ = argmax(1 .- pairwise(cosine_dist, SPopObj', W), dims = 2)
            Region_S= map(c -> c.I[2], Region_S_[:,1])

            #Sort subregions ascend
            Regions = sort(Region_S)
            flag = maximum(Regions)

            #Counters
            counters = [ count(r -> r == i, Regions) for i in 1:flag]
            # countmap(Region_S)
            
            #Count the subregions with more indiviudals
            #most_crowded=filter(r -> r == flag, Region_S)
            most_crowded = argmax(counters)
                        
            #S_crowdest is the set of individuals from the most crowded subregion
            S_crowdest = S[Regions .== most_crowded] # map(x -> x.rank .== Region_S[most_crowded]) 
            #S_crowdest=S(Region==most_crowded);                 % S_crowdest is the set of individuals from the most crowded subregion

            dist = pairwise(Euclidean(), fvals(S_crowdest)',fvals(S_crowdest)', dims=2)[:,1]
            dist[dist .== 0] .= Inf
            # dist = map(x->  x == 0 ? Inf : x , dist)

            # r = find(min(min(dist))==dist);
            r = findall(dist .== minimum(dist));
            
            #St is the set of individuals having the smallest distance in S_crowdest
            St=S_crowdest[r];   
          
            Region_St_ = argmax(1 .- pairwise(cosine_dist, fvals(St)',W), dims = 2)
            Region_St = map(c -> c.I[2], Region_St_[:,1])
            
            Z = ideal(fvals(St))
    
            g_tch = maximum(abs.(fvals(St) .- zideal') ./ W[:, Region_St]', dims = 2)[:,1]
            order = argmax(g_tch)
            x_worst=St[order]
            i_worst = findfirst(s -> s === x_worst, S)
            deleteat!(S, i_worst)
            #S = setdiff(S, [x_worst]) 
            # S = S[S .∉ Ref(x_worst)]   
        end
        CA = S
        return CA


    #3rd condition
    else

        #Obtain non feasible solutions
        # SI = Hc[Hc .∉ Ref(Sc)]
        SI = setdiff(Hc, Sc)

        # differs from platemo
        NSIObj = fvals(SI)
        NSIObj = ( NSIObj .- ideal(NSIObj)') ./ (nadir(NSIObj)' - ideal(NSIObj)')

        f1 = Metaheuristics.sum_violations.(SI) 
        
        # Tchevychev decomposition
        Region_SI_ = argmax(1 .- pairwise(cosine_dist, NSIObj', W), dims = 2)
        Region_SI = map(c -> c.I[2], Region_SI_[:,1])
        
        #Ideal point
        zideal = ideal(fvals(SI))  
        f2 = maximum(abs.(fvals(SI) .- zideal') ./ W[:, Region_SI]', dims = 2)[:,1]

        # fix nans
        f2[findall(isnan, f2)] .= 1 / eps()
        
        
        #Non dominated sorting of non feasible solutions

        #Associate each objective function to each nonfeasible solution

        #Size of nonfeasible (|si| = |f1| = |f2|)
        nonf = length(SI)
        for x in 1:nonf
            SI[x].f[1] = f1[x]
            SI[x].f[2] = f2[x]
        end

        #Nondominated ranking with modified solutions
        Metaheuristics.fast_non_dominated_sort!(SI)
        FrontNO = map(s -> s.rank, SI)
        
        #Join Sc with S
        S = vcat(Sc,S)

        #|Sc| < N, Sc U Fi ...N
        last = 0
        MaxFronts = maximum(FrontNO)
        for i in 1:MaxFronts
            fpj = filter(s -> s.rank == FrontNO[i],SI)
            append!(S, fpj)

            if length(S) >= N
                last = i
                break
            end
        end
        
        #|S| > N
        fp = filter(s -> s.rank == last,SI)
        delete_n = length(S) - N;
   
        #Delete solutions
        while delete_n > 0
            cvs_max = maximum(x->x.sum_violations, fp);
            xw = filter(s -> s.sum_violations == cvs_max,fp)
            S = S[S .∉ Ref(xw)]  
            delete_n -= 1;  
        end

        return S
    end

    return CA
end
