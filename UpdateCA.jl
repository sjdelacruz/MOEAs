using Distances
using Metaheuristics


function UpdateCA(CA,Q,W)

    # Initial structures
    
    #Output
    S=[] 

    #Feasible solutions
    Sc = [] 

    #Number of solutions (at the beginning is the initial population)
    N=size(Q,1)


    i=1
    Hc = vcat(CA,Q)


    #Feasible solutions 
    Sc = filter((x) -> x.is_feasible == true, Hc)

    #Number of feasible solutions   
    feasibles = size(Sc,1)
    
    if feasibles == N
        CA = Sc
    elseif feasibles >= N
        
        #Non dominated sorting to Sc
        Metaheuristics.fast_non_dominated_sort!(Sc)
        FrontNO = map(s -> s.rank, Sc)

        #Size of S
        zpoints = size(S,1)

        while zpoints < N
            S = append!(S, Si[FrontNO[i].fvals])
            i = i+1
            zpoints = size(S,1)
        end

        if zpoints > N

            #Ideal and nadir points
            zideal = ideal(S) 
            znad = nadir(S)
            
            #Normalize_points
            normalize = map((x) ->  (x - zideal)/(znad-zideal), S)

            
            # association of S and W 
            Region_CA_ = argmax(1 .- pairwise(cosine_dist, fvals(CA)', W), dims = 2)
            Region_CA = map(c -> c.I[2], Region_CA_[:,1])

        end
        CA = S;


    else

        #Obtain elements 
        Si = Hc[Hc .∉ Ref(Sc)]

        #Non dominated sorting
        Metaheuristics.fast_non_dominated_sort!(Si)
        FrontNO = map(s -> s.rank, Si)

        #|Sc| < N, Sc U Fi ...N
        
        while feasibles < N
            append!(S, Si[FrontNO[i]])
            i = i+1
        end

        #Size of S
        zpoints = size(S,1)
        
        #Number of fronts
        last = size(FrontNO)
        
        #|S| > N
        while zpoints > N
            xw = argmax(Si[FrontNO[i-1]])
            S = S[S .∉ Ref(xw)]    
        end
        CA = S
    end

    """

        

        #Fill feasible solutions
        fcounter = size(Sc,1)
        i=0
        while fcounter < N
            
            S = vcat(Sc,FrontNO[i])
            fcounter += size(FrontNo[i],1)
            i=i+1
        end

        if size(S,1) > N

            zideal = ideal(Hc)
            znad = nadir(Hc)

            #Objective values normalization
            for x in S
                x.objs = (x.objs - zideal) / znad - zideal
            end 
        
        
            #Asociation
            Region_Hc_ = argmax(1 .- pairwise(cosine_dist, fvals(Hc)', W), dims = 2)
            Region_Hc = map(c -> c.I[2], Region_Hc_[:,1])
            
            #Associat the individuals in Hd and CA with their corresponding subregions
            Region_CA_ = argmax(1 .- pairwise(cosine_dist, fvals(CA)', W), dims = 2)
            Region_CA = map(c -> c.I[2], Region_CA_[:,1])

            #Find the most crowded subregion 
            crowded_region = [];
            dist_ind = []
            
            #Adjust distances
            #for x in crowded_region
            #end

            #Minimun values
            St = argmin(dist)
        end
        
        Si = 
        CA = S
    else feasibles > N

        
    end
    """
end
