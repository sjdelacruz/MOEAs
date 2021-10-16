using Distances
using Metaheuristics


function UpdateDA(CA,Q,W,N)

    # Initial structures
    S=[]
    Sc = []
    i=1
    Hc = vcat(CA,Q)
    N=size(W,1)

    #Feasible solutions
    Sc = map(x ->  Metaheuristics.violationsSum(x) == 0, Hc)

    #Number of feasible solutions
    feasibles = size(Sx,1)

    if feasibles == N
        CA = Sc
    elseif 

        #Non dominated sorting
        Metaheuristics.fast_non_dominated_sort!(Sc)
        FrontNO = map(s -> s.rank, Sc)

        #Fill feasible solutions
        fcounter = size(Sc,1)
        i=0
        while fcounter < N
            S = vcat(Sc,FrontNO[i])

            fcounter += size(FrontNo[i],1)
            i++
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
            for x in crowded_region
            end

            #Minimun values
            St = argmin(dist)
        end
        
        Si = 
        CA = S
    else feasibles > N

        #Obtain elements 
        Si = Hc[(!in).(Hc,Ref(Sc))]

        #Non dominated sorting
        Metaheuristics.fast_non_dominated_sort!(Si)
        FrontNO = map(s -> s.rank, Si)

        ndif = size(Sc,1)
        while ndif < N
            S = vcat(S,FrontNO[i])
            i++
            ndif += size(FrontNO[i],1)
        end

        elems = size(S,1)

        while elems > N
            xw = argmax()
        end

    end


    