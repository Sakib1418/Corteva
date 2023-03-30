def [Population,Params] = Main(varargin)
format long g

## Call initialization defs

#Variable input
#varargin can contain the structure array Params, if it has been previously
#created. This saves considerable time from having to read the input files
#and create a new structure array. 

if ~isempty(varargin)
    Params = varargin{1}
else
    #Get user defined parameters and store them to a temporary variable to be
    #replaced after the prepareParams() def.
    tmp = readData()

    #Call the prepareParams() def to create a permanent params array and
    #initilize the population variables.
    Params = prepareParams_2(tmp)


## Initialize life-cycle storage variables

#Create a structure array to hold the various populations
Population = struct

#Create (redundant) frequently used variables to make code more readable
nCohorts = Params.General.nCohorts         #number of cohorts to simulate
nYears = Params.General.nYears             #number of years to simulate
nLoci = Params.General.nLoci               #number of gene loci under consideration
nGeno = Params.General.nGeno               #number of unique genotypes under consideration
fieldSize = Params.General.fieldSize       #field are in square meter

#Structure array for storing various life-cycle stages
#seed bank
Population.seedBank = np.zeros(nYears,nGeno,nCohorts)

#Lower seed bank density
Population.lowerBank = np.zeros(nYears,nGeno)

#immigrant seed
Population.immSeed = np.zeros(nYears,nGeno,nCohorts)

#emigrant seed
Population.emSeed = np.zeros(nYears,nGeno,nCohorts)

#germination storage
Population.germination = np.zeros(nYears,nGeno,nCohorts)

#ungerminated storage
Population.ungerminated = np.zeros(nYears,nGeno,nCohorts)

#herbicide storage
Population.herbicide1 = np.zeros(nYears,nGeno,nCohorts)
Population.herbicide2 = np.zeros(nYears,nGeno,nCohorts)
Population.herbicide3 = np.zeros(nYears,nGeno,nCohorts)
Population.herbicide4 = np.zeros(nYears,nGeno,nCohorts)

#adult plant storage
Population.mature = np.zeros(nYears,nGeno,nCohorts)

#hand removal storage
Population.hand = np.zeros(nYears,nGeno,nCohorts)

#immigrant pollen
Population.immPol = np.zeros(nYears,nGeno,nCohorts)

#emigrant seed
Population.emPol = np.zeros(nYears,nGeno,nCohorts)

#seed produced storage
Population.newSeed = np.zeros(nYears,nGeno,nCohorts)

#seed population after mutation
Population.mutated = np.zeros(nYears,nGeno,nCohorts)

#seed population after predation
Population.predation = np.zeros(nYears,nGeno,nCohorts)

#seed population after winter (note this happens after all cohorts
#have done their thing, so cohort data isn't stored)
Population.winter = np.zeros(nYears,nGeno)

#Resistance metrics data
Population.Res = np.zeros(1, 3)

## Final initialization steps

#Put the SeedBank array from the initialization defs into the
#population.seedBank array for year 1, all genotypes, and cohort 1
Population.seedBank(1,:,1) = Params.General.UpperSeedBank
Population.lowerBank(1,:) = Params.General.LowerSeedBank


## MAIN SIMULATION

#Run simulation of nCohorts and nYears, storing information at each step in
#the simulation NOTE: (1 x nGeno) variables called curPop (for current
#population) and curBank (for current seed bank) will be used to store the
#population in the current loop

#Each stage in the life cycle is given its own section, denoted by ## and a
#horizontal line if editing in Matlab. The format of the code for each
#life-cycle stage is:
#A. Get parameters relevant for the current stage 
#B. Call transition def for the current stage
#C. Adjust populations if needed (e.g. extinction)
#D. Store current population in population structure array

for y = 1:nYears
    tempSeedBank = 0

    #If it till frequency = current year, upper and lower seed banks
    #will be swapped, otherwise, nothing changes
    [Population.seedBank(y,:,1),Population.lowerBank(y,:)] = ...
        deepTill(Population.seedBank(y,:,1),Population.lowerBank(y,:),Params.General.TillingFreq,y)

    for c = 1:nCohorts
        ## Curent seed bank + immigration

        #Allow seed immigration from neighboring field
        ImmSeed =  seedImmigration(nGeno)

        #Set curPop = to the seedBank + immigration
        seedBank = Population.seedBank(y,:,c) + ImmSeed

        #Call extinction def 
        seedBank = extinction(seedBank,fieldSize)

        ## Germination (starts with SeedBank)

        #A. Get parameters

        #germination fraction for current cohort/year
        fGerm = Params.Germination.survivalFraction(y,c)

        #B. Call transition def

        #survival def with germination parameters (storing
        #population in a variable curPop (for current population)
        germination = survival(seedBank,fGerm)

        #C. Adjust populations 

        #subtract germinated seed from the SeedBank and stored as
        #ungerminated. Ungerminated seed gets added to next cohort's seed
        #bank (either this year or next year)
        ungerminated = seedBank - germination

        #store curPop value in population structure
        germination = extinction(germination,fieldSize)

        ## Cultivation (survival def)

        #survival fraction for current cohort/year (fraction surviving
        #predation)
        fCult = Params.Cultivation.survivalFraction(y,c)

        #call survival def for predation
        cult = survival(germination,fCult)

        #store seed remaining after predation
        cult = extinction(cult,fieldSize)

        ## Herbicide 1 application

        #survival fraction for current cohort/year for each genotype
        fSS = 1 - Params.Herbicide1.efficacySS(y,c)
        fRS = 1 - Params.Herbicide1.efficacyRS(y,c)
        fRR = 1 - Params.Herbicide1.efficacyRR(y,c)

        #target locus
        targetLocus = Params.Herbicide1.targetLocus

        #call selective survival def for herbicide1
        herbicide1 = selectiveSurvival(cult,[fSS,fRS,fRR],targetLocus,Params.key)

        #store value in population structure 
        herbicide1 = extinction(herbicide1,fieldSize)

        ## Herbicide 2 application

        #survival fraction for current cohort/year for each genotype
        fSS = 1 - Params.Herbicide2.efficacySS(y,c)
        fRS = 1 - Params.Herbicide2.efficacyRS(y,c)
        fRR = 1 - Params.Herbicide2.efficacyRR(y,c)

        #target locus
        targetLocus = Params.Herbicide2.targetLocus

        #call selective survival def for herbicide2
        herbicide2 = selectiveSurvival(herbicide1,[fSS,fRS,fRR],targetLocus,Params.key)

        #store value in population structure 
        herbicide2 = extinction(herbicide2,fieldSize)

        ## Herbicide 3 application

        #survival fraction for current cohort/year for each genotype
        fSS = 1 - Params.Herbicide3.efficacySS(y,c)
        fRS = 1 - Params.Herbicide3.efficacyRS(y,c)
        fRR = 1 - Params.Herbicide3.efficacyRR(y,c)

        #target locus
        targetLocus = Params.Herbicide3.targetLocus

        #call selective survival def for herbicide3
        herbicide3 = selectiveSurvival(herbicide2,[fSS,fRS,fRR],targetLocus,Params.key)

        #store value in population structure 
        herbicide3 = extinction(herbicide3,fieldSize)

        ## Herbicide 4 application

        #survival fraction for current cohort/year for each genotype
        fSS = 1 - Params.Herbicide4.efficacySS(y,c)
        fRS = 1 - Params.Herbicide4.efficacyRS(y,c)
        fRR = 1 - Params.Herbicide4.efficacyRR(y,c)

        #target locus
        targetLocus = Params.Herbicide4.targetLocus

        #call selective survival def for herbicide4
        herbicide4 = selectiveSurvival(herbicide3,[fSS,fRS,fRR],targetLocus,Params.key)

        #store value in population structure 
        herbicide4 = extinction(herbicide4,fieldSize)

        ## Competition to maturity

        #survival fraction for current cohort/year for each genotype
        A = Params.Mature.maxPlants(y,c)
        B = Params.Mature.cropParam(y,c)
        C = Params.Mature.weedCompetition(y,c)

        #call competition def for herbicide4
        mature = competition(herbicide4,A,B,C)

        #store curPop value in population structure 
        mature = extinction(mature,fieldSize)

        ## Hand Removal (survival def)

        #survival fraction for current cohort/year (fraction surviving
        #hand removal)
        fHand = Params.Hand.survivalFraction(y,c)

        #call survival def for predation
        hand = survival(mature,fHand)

        #store seed remaining after predation
        hand = extinction(hand,fieldSize)

        ## Seed production competition and Mating

        #survival fraction for current cohort/year for each genotype
        A = Params.SeedProd.maxYield(y,c)
        B = Params.SeedProd.cropParam(y,c)
        C = Params.SeedProd.weedCompetition(y,c)

        #call selective survival def for herbicide1
        seedYield = competition(sum(hand),A,B,C)

        # Call mating def
        newFrac = mating(hand,Params.General.selfingCoeff,Params.key)

        # Newly produced seed
        seedProd = newFrac * seedYield * Params.General.femaleFrac
        seedProd = extinction(seedProd,fieldSize)

        ## Mutation
        #Called for loci 1-4
        #Has a bypass built in to avoid repetitive def calls if
        #mutation parameters are all zero
        
        #Locus 1
        #mutation parameters for current cohort/year
        targetLocus = 1
        
        muS2R = Params.Mutation.muS2R_locus1(y,c)
        muR2S = Params.Mutation.muR2S_locus1(y,c)

        if muS2R==0 && muR2S==0
            mutated_locus1 = seedProd
        else
            #mutation def
            mutated_locus1 = mutation(seedProd,muS2R,muR2S,targetLocus,Params.key)
            mutated_locus1 = extinction(mutated_locus1,fieldSize)
        
        
        #Locus 2
        #mutation parameters for current cohort/year
        targetLocus = 2
        muS2R = Params.Mutation.muS2R_locus2(y,c)
        muR2S = Params.Mutation.muR2S_locus2(y,c)

        if muS2R==0 && muR2S==0
            mutated_locus2 = mutated_locus1
        else
            #mutation def
            mutated_locus2 = mutation(mutated_locus1,muS2R,muR2S,targetLocus,Params.key)
            mutated_locus2 = extinction(mutated_locus2,fieldSize)
        
        
        #Locus 3
        #mutation parameters for current cohort/year
        targetLocus = 3
        muS2R = Params.Mutation.muS2R_locus3(y,c)
        muR2S = Params.Mutation.muR2S_locus3(y,c)

        if muS2R==0 && muR2S==0
            mutated_locus3 = mutated_locus2
        else
            #mutation def
            mutated_locus3 = mutation(mutated_locus2,muS2R,muR2S,targetLocus,Params.key)
            mutated_locus3 = extinction(mutated_locus3,fieldSize)
        
        
        #Locus 4
        #mutation parameters for current cohort/year
        targetLocus = 4
        muS2R = Params.Mutation.muS2R_locus4(y,c)
        muR2S = Params.Mutation.muR2S_locus4(y,c)

        if muS2R==0 && muR2S==0
            mutated_locus4 = mutated_locus3
        else
            #mutation def
            mutated_locus4 = mutation(mutated_locus3,muS2R,muR2S,targetLocus,Params.key)
            mutated_locus4 = extinction(mutated_locus4,fieldSize)
        

        ## Predation (survival def)

        #survival fraction for current cohort/year (fraction surviving
        #predation)
        fPred = Params.Predation.survivalFraction(y,c)

        #call survival def for predation
        predation = survival(mutated_locus4,fPred)

        #store seed remaining after predation
        predation = extinction(predation,fieldSize)

        ## Sort new seed
        #Determine where the newly produced seed goes. Two options
        if Params.General.seedDelay == 0
            #Put newly produced seed in next cohort's seed bank
            if c < nCohorts
                Population.seedBank(y,:,c+1) = Population.seedBank(y,:,c+1) + predation + ungerminated
            else
                tempSeedBank = tempSeedBank + ungerminated + predation
            

        else
            #Put newly produced seed in temporary seed bank, to be added to
            #next year's seed bank.
            tempSeedBank = tempSeedBank + predation
            if c < nCohorts
                Population.seedBank(y,:,c+1) = Population.seedBank(y,:,c+1) + ungerminated
            else
                tempSeedBank = tempSeedBank + ungerminated
            
        

        ## Store information from this iteration of the loop in population array
        Population.immSeed(y,:,c) = ImmSeed
        Population.ungerminated(y,:,c) = ungerminated
        Population.germination(y,:,c) = germination
        Population.herbicide1(y, :, c) = herbicide1
        Population.herbicide2(y, :, c) = herbicide2
        Population.herbicide3(y, :, c) = herbicide3
        Population.herbicide4(y, :, c) = herbicide4
        Population.mature(y,:,c) = mature
        Population.hand(y,:,c) = hand
        Population.seedProd(y,:,c) = seedProd
        Population.predation(y,:,c) = predation

    
# 
#     #Update weed resistance data according to this year's population.
#     Population.Res = resMetrics(Population.Res, y, Population.seedBank, Params.key, ...
#                                 Params.General.latentTime, Params.General.criticalResTime, ...
#                                 Params.General.completeResTime)

    ## Winter (survival def)
    #If in the last year, break out of loop
    if y == nYears
        break
    else
        #survival fraction for current cohort/year (fraction surviving
        #predation)
        fWint_upper = Params.Winter.uppersurvivalFraction(y)
        fWint_lower = Params.Winter.lowersurvivalFraction(y)

        #call winter survival def for upper seed bank
        winter = survival(tempSeedBank,fWint_upper)
        winter = extinction(winter,fieldSize)
        
        Population.winter(y,:) = winter

        #store current population in next years seed bank     
        Population.seedBank(y+1,:,1) = Population.seedBank(y+1,:,1) + winter

        #lower seed bank
        Population.lowerBank(y+1,:) = survival(Population.lowerBank(y,:),fWint_lower)
        Population.lowerBank(y+1,:) = extinction(Population.lowerBank(y+1,:),fieldSize)
    


plotPopulation(Population.seedBank,Params.key,'total')


