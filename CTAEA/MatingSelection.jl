using Metaheuristics

function MatingSelection(S)
    # Mating selection of C-TAEA

    number=length(S)
    rnd = rand(1:number,2)
    x_1=rnd[1]
    x_2=rnd[2]
    CV1 = Metaheuristics.sum_violations(S[x_1]) # sum(max(0,S(x_1).con),2);
    CV2 = Metaheuristics.sum_violations(S[x_2]) # sum(max(0,S(x_2).con),2);

    if CV1 > 0 && CV2 > 0
        return S[rand(1:number)]
    end

    # tournament
    if Metaheuristics.compare(S[x_1], S[x_2]) == 1
        return S[x_1]
    end

    return S[x_2]
end

#=
function test()
    _, _, population = Metaheuristics.TestProblems.get_problem(:ZDT3)

    MatingSelection(population)
end

test()
=#
