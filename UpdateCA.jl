using Distances
using Metaheuristics

#N Number of solutions (at the beginning is the initial population)

function UpdateCA(CA,Q,W,N)
    
    #Output
    S=empty(Q)
    i=1

    #Hybrid population
    Hc = vcat(CA,Q)
     
    #Feasible solutions 
    Sc = filter((x) -> x.is_feasible == true, Hc)

    #Number of feasible solutions   
    feasibles = size(Sc,1)
    
    #1st condition
    if feasibles == N
        CA = Sc

    # 2nd condition
    elseif feasibles >= N
        
        #Non dominated sorting to Sc
        Metaheuristics.fast_non_dominated_sort!(Sc)
        FrontNO = map(s -> s.rank, Sc)

        for i in  nfronts
            fp = findall(s -> s.rank == FrontNO[i],Sc)
            S = append!(S, fp)

            zs = size(S,1)
            if zs > N
                break
            end
        end

        #Size of S
        zs = size(S,1)

        while zs > N

            #Ideal and nadir points
            zideal = ideal(S) 
            znad = nadir(S)
            
            #Normalize_points
            map((x) ->  x.f = (x.f .- zideal) / (znad .-zideal), S)

            
            # association of S and W 
            Region_S_ = argmax(1 .- pairwise(cosine_dist, S, W), dims = 2)
            Region_S= map(c -> c.I[2], Region_S_[:,1])

            #Sort subregions ascend
            Regions = sort(Region_S)

            #Counters
            counters = map(x -> x.elems, Region_S)

            #Count the subregions with more indiviudals
            most_crowded=max(counter);
            
            #S_crowdest is the set of individuals from the most crowded subregion
            S_crowdest=map(x -> x.rank == most_crowded) 

            dist = pairwise(Euclidean, vals(S_crowdest),fvals(S_crowdest));                     
            dist = map(x->  x == 0 ? Inf16 : 0 , dist)

            r = find(min(min(dist))==dist);
            
            #St is the set of individuals having the smallest distance in S_crowdest
            St=S_crowdest(r);   
          
            Region_St = argmax(1 .- pairwise(cosine_dist, fvals(St),W), dims = 2)

            Z = ideal(fvals(St))
    
            g_tch=maximum(abs.(fvals(St) .- zideal) ./ W[:, Region_St], dims = 2)[:,1]
            order=argmax(g_tch)
            x_worst=St[order]
            S = S[S .∉ Ref(x_worst)]
        
        end

        CA = S;

    #3rd condition
    else

        #Obtain non feasible solutions
        SI = Hc[Hc .∉ Ref(Sc)]
        f1= map(x -> max(0, x), Metaheuristics.sum_violations.(SI));

        Region_SI_ = argmax(1 .- pairwise(cosine_dist, fvals(SI),W), dims = 2)
        Region_SI = map(c -> c.I[2], Region_SI_[:,1])
    
        #Ideal point
        zideal = ideal(fvals(SI))
        subp = map(x -> broadcast(abs, x - zideal), fvals(SI))
        
        f2 = argmax(subp ./ W[Region_SI,:], dims = 2);
        
        #Non dominated sorting of non feasible solutions
        Metaheuristics.fast_non_dominated_sort!(SI)
        FrontNO = map(s -> s.rank, SI)
        MaxFronts = size(FrontNO,1)

        # join Sc with S
        S = vcat(Sc,S)

        #|Sc| < N, Sc U Fi ...N

        last = 0
        for i in MaxFronts
            fp = findall(s -> s.rank == FrontNO[i],SI)
            append!(S, fp)
            i = i+1 

            nfpoints = size(S,1)
            if nfpoints > N
                last = size(FrontNO-1)
                break
            end
        end
        

        #|S| > N
        fp = findall(s -> s.rank == last,SI)
        delete_n=size(S,2)-N;
   
        while delete_n > N

            xw = maximum(x->x.sum_violations, fp)
            S = S[S .∉ Ref(xw)]  
            delete_n = delete_n-1;  
        end

        CA = S
    end

    return CA
end

function count_by_group(regions)

    Counts = []

    for i in regions
        
        x = counter(i > i == 1, )
    end
end