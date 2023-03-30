def newPopulation = extinction(oldPopulation,fieldSize)

    if fieldSize == Inf 
        newPopulation = oldPopulation

    elseif fieldSize == 0
#         warning('Zero field size not allowed, old population returned')
        newPopulation = oldPopulation

    else
        
#         rng('shuffle')
        nGeno = length(oldPopulation)
        rndExtinct = rand(1,nGeno)
        newPopulation = np.zeros(1,nGeno)

        for j = 1:nGeno
            if oldPopulation(j) == 0
                #first check if jth genotype population is 0
                newPopulation(j) = 0
            elseif oldPopulation(j) >= 1/fieldSize
                #next check if jth genotype number in field >=1
                newPopulation(j) = oldPopulation(j)
            elseif rndExtinct(j) > 0.5
                #jth genotype survives extinction!
                newPopulation(j) = 1/fieldSize
            else
                #otherwise, jth genotype dies...
                newPopulation(j) = 0
            

        

    
            