using Package
using Distances
using Metaheuristics
using StatsBase

#N Number of solutions (at the beginning is the initial population)

function UpdateCA(CA,Q,W_,N)
    
    #Change to dimensions weighted vectors
    W = hcat(W_...)

    #Output
    S=empty(Q)

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

        nfront = length(FrontNO)

        #Fill S with fronts i=1 ... nfront
        for i in  range(1, length= nfront)
            fp = filter(s -> s.rank == FrontNO[i],Sc)
            S = append!(S, fp)
            if size(S,1) >= N
                break
            end
        end

        while size(S,1) > N

            #Ideal and nadir points
            zideal = ideal(fvals(S))' 
            znad = nadir(fvals(S))'
            
            #Normalize_points
            map(x ->  x.f = ((x.f .- zideal) ./ (znad .-zideal)), S)

            # association of S and W 
            Region_S_ = argmax(1 .- pairwise(cosine_dist, fvals(S)', W), dims = 2)
            Region_S= map(c -> c.I[2], Region_S_[:,1])

            #Sort subregions ascend
            Regions = sort(Region_S)

            #Counters
            counters = countmap(Region_S)
            
            flag = max(counters)
            #Count the subregions with more indiviudals
            most_crowded=filter(r -> r = flag,Region_S)
            
            #S_crowdest is the set of individuals from the most crowded subregion
            S_crowdest = map(x -> x.rank == Region_S[most_crowded]) 

            dist = pairwise(Euclidean, fvals(S_crowdest)',fvals(S_crowdest)');                     
            dist = map(x->  x == 0 ? Inf16 : 0 , dist)

            r = find(min(min(dist))==dist);
            
            #St is the set of individuals having the smallest distance in S_crowdest
            St=S_crowdest[r];   
          
            Region_St_ = argmax(1 .- pairwise(cosine_dist, fvals(St)',W), dims = 2)
            Region_St = map(c -> c.I[2], Region_St_[:,1])
            
            Z = ideal(fvals(St))'
    
            g_tch=maximum(abs.(fvals(St)' .- zideal) ./ W[:, Region_St], dims = 2)[:,1]
            order=argmax(g_tch)
            x_worst=St[order]
            S = S[S .∉ Ref(x_worst)]   
        end
        CA = S;

    #3rd condition
    else

        #Obtain non feasible solutions
        SI = Hc[Hc .∉ Ref(Sc)]

        f1= Metaheuristics.sum_violations.(SI)
        
        Region_SI_ = argmax(1 .- pairwise(cosine_dist, fvals(SI)', W), dims = 2)
        Region_SI = map(c -> c.I[2], Region_SI_[:,1])
        
        #Ideal point
        zideal = ideal(fvals(SI))'
        
        f2 = maximum(abs.(fvals(SI) .- zideal) ./ W[:, Region_SI]', dims = 2)[:,1]
        
        #Non dominated sorting of non feasible solutions

        #Associate each objective function to each nonfeasible solution

        #Size of nonfeasible (|si| = |f1| = |f2|)
        nonf = length(SI)
        for x in range(1,length=nonf)
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
        MaxFronts = length(FrontNO)
        for i in range(1, length = MaxFronts)
            fpj = filter(s -> s.rank == FrontNO[i],SI)
            append!(S, fpj)

            if size(S,1) >= N
                last = i
                break
            end
        end
        
        #|S| > N
        fp = filter(s -> s.rank == last,SI)
        delete_n=size(S,1)-N;
   
        #Delete solutions
        while delete_n > 0
            cvs_max = maximum(x->x.sum_violations, fp);
            xw = filter(s -> s.sum_violations == cvs_max,fp)
            S = S[S .∉ Ref(xw)]  
            delete_n = delete_n-1;  
        end

        CA = S
    end

    return CA
end
