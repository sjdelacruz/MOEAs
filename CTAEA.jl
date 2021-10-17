using Metaheuristics
##Duda al llamar la funcion Algorithm

mutable struct CTAEA <: Metaheuristics.AbstractParameters
    nobjectives::Int
    N::Int #individuos
end

function CTAEA(weights; N = 100, information = Information(), options = Options())
    
    if isempty(weights)
        error("Provide weighted vectors")
    end

    nobjectives = length(weights[1])

    parameters = CTAEA(nobjectives,N)

    alg = Algorithm(
        parameters,
        information = information,
        options = options,)
    
    alg
end

function initialize!(
    status,
    parameters::CTAEA,
    problem::Metaheuristics.AbstractProblem,
    information::Information,
    options::Options,
    args...;
    kargs...
)
        
    if options.iterations == 0
        options.iterations = 500
    end

    if options.f_calls_limit == 0
        options.f_calls_limit = options.iterations * parameters.N + 1
    end

    status = gen_initial_state(problem,parameters,information,options,status)
    D = size(problem.bounds, 2)
    parameters.nobjectives = length(status.population[1].f)

    return status
end

function update_state!(
    status::State,
    parameters::CTAEA,
    problem::Metaheuristics.AbstractProblem,
    information::Information,
    options::Options,
    args...;
    kargs...)

    N = parameters.N
    D = size(problem.bounds, 2)

end
    
    """CA=UpdateCA([],status.population,parameters.Wv);            
    DA=UpdateDA(CA,[],status.population,parameters.Wv);
    
    
    stop_criteria!(status, parameters, problem, information, options)
    while !status.stop    
    
        # mating pool choosing
        #calculate the ratio of non-dominated solutions of CA and DA in Hm
        Hm=[CA,DA];                         
        [FrontNo,~]=NDSort(Hm.objs,inf);
        FrontNo_C=FrontNo(1:ceil(length(Hm)/2));
        Nc=size(find(FrontNo_C==1),2);      
        Pc=Nc/length(Hm);
        FrontNo_D=FrontNo(ceil(length(Hm)/2)+1:length(Hm));
        Nd=size(find(FrontNo_D==1),2);      
        Pd=Nd/length(Hm);

        # calculate the proportion of non-dominated solutions in CA
        [FrontNo,~]=NDSort(CA.objs,inf);
        NC=size(find(FrontNo==1),2); 

        # PC denotes the proportion of non-dominated solutions in CA,it is different from Pc     
        PC=NC/length(CA);                     

        #reproduction
        Q=[];
        for i=1:size(Wv,1)
            if Pc>Pd
                P1=MatingSelection(CA); 
            else
                P1=MatingSelection(DA);
            end
            pf=rand();

            if pf<PC
                P2=MatingSelection(CA);
            else
                P2=MatingSelection(DA);
            end

            MatingPool=[P1,P2];
            
            #SBX and polynomial mutation
            Offspring=OperatorGAhalf(MatingPool);

            Q=[Q,Offspring];
        end

    # update CA and DA
        CA=UpdateCA(CA,Q,W);
        DA=UpdateDA(CA,DA,Q,W);
    end
    """

function stop_criteria_ctaea(
    status::State,
    parameters::CTAEA,
    problem::Metaheuristics.AbstractProblem,
    information::Information,
    options::Options,
    args...;
    kargs...
    )

    return status.iteration > options.iterations
end

function final_stage!(
    status,
    parameters::CTAEA,
    problem::Metaheuristics.AbstractProblem,
    information::Information,
    options::Options,
    args...;
    kargs...
)
    status.final_time = time()
end


#######################################################################################################################################
"""
function UpdateCA(CA, Q, Wv)
    
    S=[];       # S is the set used for output
    Sc=[];      # Sc is used to collect feasible solutions
    Hc=[CA,Q];   #Feasibility population and offspring
    N=size(Wv,1);

    p1=sum(max(0,Hc.g),2);   
    p2=sum(max(0,Hc.h),2);  
    CV= p1+p2;  

    Sc=[Sc,Hc(CV==0)];

    if length(Sc)==N
        UpdatedCA=Sc;
    elseif length(Sc)>N

        [FrontNO,MaxNO]=get_non_dominated_solutions(population) ##Verificar metodo
        for i=1:MaxNO
            S=cat(2,S,Sc(FrontNO==i));
            if length(S)>=N
                break;
            end
        end

        while length(S)>N 
            
            #normalization
            Zmax=max(S.f,[],1);
            Zmin=min(S.f,[],1);

            SPopObj=(S.f-repeat(Zmin,size(S.f,1),1))./(repmat(Zmax,size(S.f,1),1)-repmat(Zmin,size(S.f,1),1));

            # associate each solution in S with their corresponding subregion
            [~,Region] = max(1-Euclidean(SPopObj,Wv),[],2)
            [value,~]=sort(Region,'ascend');
            flag=max(value);

            #counter denotes the number of indiviudals in each subregion
            counter=histc(value,1:flag);                        
            [~,most_crowded]=max(counter);

            # S_crowdest is the set of individuals from the most crowded subregion
            S_crowdest=S(Region==most_crowded);                 
            dist=pdist2(S_crowdest.objs,S_crowdest.objs);                     
            dist(dist==0)=inf;
            [row,~]=find(min(min(dist))==dist);

            # St is the set of individuals having the smallest distance in S_crowdest
            St=S_crowdest(row);                                 
            [~,Region_St] = max(1-pdist2(St.objs,W,'cosine'),[],2);
            Z = min(St.objs,[],1);
            g_tch=max(abs(St.objs-repmat(Z,length(St),1))./W(Region_St,:),[],2);
            [~,order]=max(g_tch);
            x_wrost=St(order);
            S=setdiff(S,x_wrost);
        end
        UpdatedCA=S;

    elseif length(Sc)<N
        # SI is the set of infeasible solutions in Hc
        SI=setdiff(Hc,Sc);	
        f1=sum(max(0,SI.cons),2);
        [~,Region_SI] = max(1-pdist2(SI.objs,W,'cosine'),[],2);
        Z = min(SI.objs,[],1) ;
        f2=max(abs(SI.objs-repmat(Z,length(SI),1))./W(Region_SI,:),[],2);
        PopObj=[f1,f2];           
        [FrontNO,MaxNO]=NDSort(PopObj,inf);                
        S=[S,Sc];
        for i=1:MaxNO
            S=cat(2,S,SI(FrontNO==i));
            if length(S)>=N 
                last=i;
                break;
            end
        end

        # find the individuals in the last front joined into S
        F_last=SI(FrontNO==last);	
        delete_n=size(S,2)-N;
        CV=sum(max(0,F_last.cons),2);
        [~,index]=sort(CV,'descend');
        x_wrost=F_last(index(1:delete_n));
        S=setdiff(S,x_wrost);
        UpdatedCA=S;
    end

end

"""

