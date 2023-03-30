def mutatedPopulation = mutation(population,muS2R,muR2S,targetLocus,key)

mu = muR2S
nu = muS2R

#define a new population array and initilize it to the current population
mutatedPopulation = population

#cycle through genotypes
for j = 1:size(key,2)
    #Note: Since mutation is taken one locus at a time, each
    #genotype can only become 2 other genotypes and vice versa, i.e.
    #for locus 1, [1,2,1] -> [0,2,1] or [2,2,1] only. For each
    #genotype, we need to first locate the other two that it can
    #interact with
    tmp = mod(j,3^(targetLocus-1))+1

    indSS = tmp
    indRS = indSS + 3^(targetLocus-1)
    indRR = indRS + 3^(targetLocus-1)

    if key(targetLocus,j) == 2
        #if the current genotype + locus has 2 resistant alleles
        mutatedPopulation(indRR) = mutatedPopulation(indRR) - 2*mu*population(indRR) + mu^2*population(indRR)
        mutatedPopulation(indRS) = mutatedPopulation(indRS) + 2*mu*population(indRR) - 2*mu^2*population(indRR)
        mutatedPopulation(indSS) = mutatedPopulation(indSS) + mu^2*population(indRR)
    elseif key(targetLocus,j) == 1
        #if the current genotype + locus has 1 resistant alleles
        mutatedPopulation(indRR) = mutatedPopulation(indRR) + nu*population(indRS)
        mutatedPopulation(indRS) = mutatedPopulation(indRS) - nu*population(indRS) - mu*population(indRS)
        mutatedPopulation(indSS) = mutatedPopulation(indSS) + mu*population(indRS)    
    else
        #if the current genotype + locus has 0 resistant alleles
        mutatedPopulation(indRR) = mutatedPopulation(indRR) + nu^2*population(indSS)
        mutatedPopulation(indRS) = mutatedPopulation(indRS) + 2*nu*population(indSS) - 2*nu^2*population(indSS)
        mutatedPopulation(indSS) = mutatedPopulation(indSS) - 2*nu*population(indSS) + nu^2*population(indSS)  
    

    population = mutatedPopulation


