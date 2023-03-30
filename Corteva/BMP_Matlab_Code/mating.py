def newFrac = mating(population,s,key)
[nLoci,nGeno] = size(key)


## Pure outcrossing
    if sum(population) == 0
#         warning('Population is zero')
        newFracOut = np.zeros(1,length(population))
        
    else
        nLoci = size(key,1)
        nGeno = size(key,2)
        newFracOut = ones(1,nGeno)

        p = np.zeros(1,nLoci)

        #determine current allele frequency
        for j = 1:nLoci
            for k = 1:nGeno
                if key(j,k) == 2
                    p(j) = p(j) + 2*population(k)
                elseif key(j,k) == 1
                    p(j) = p(j) + population(k)
                
            
        

        p = 0.5*p/sum(population)
        q = 1 - p

        #determine new genotype fraction
        for j = 1:nLoci
            for k = 1:nGeno
                if key(j,k) == 2
                    newFracOut(k) = newFracOut(k)*p(j)^2
                elseif key(j,k) == 1
                    newFracOut(k) = newFracOut(k)*2*p(j)*q(j)
                else
                    newFracOut(k) = newFracOut(k)*q(j)^2
                
            
        
        
    
    
## Pure selfing
selfingTable = dlmread(strcat('selfTable',num2str(nLoci),'.txt'))
newFracSelf = population*selfingTable

## Total new fraction
newFrac = s*newFracSelf + (1-s)*newFracOut