def p = readData( )
#Reads data from text files created by Excel and stores them into a
#structure for the program to use.

format long g

## Create struct to hold all the data.
p = struct('General',{},'Init', {},'Germination', {}, 'Cultivation', {}, ...
    'Hand', {}, 'Mature', {},'SeedProd', {},'Mutation', {}, 'Predation', {}, ...
    'Winter', {}, 'Herbicide1', {}, 'Herbicide2', {}, 'Herbicide3', {}, 'Herbicide4', {})

## Initiate static parameters.
p(1).General.currentYear = 1
p(1).General.currentCohort = 1

## Get Windows user name to construct correct file paths
username = getenv('USERNAME')

## Read and store life cycle parameter data.
path = strcat('C:\Users\', username, '\Desktop\Weed_Resistance_Data.txt')

#Read file to local variable
params_in = dlmread(path)

p(1).General.nSim = params_in(1,1)                 #number of simulations
p(1).General.nYears = params_in(2,1)               #number of years
p(1).General.nCohorts = params_in(3,1)             #number of cohorts
p(1).General.nLoci = 4                             #number of gene loci
p(1).General.nGeno = 3^p.General.nLoci             #number of genotypes
p(1).General.fieldSize = params_in(5,1)            #field size
p(1).General.selfingCoeff = params_in(6,1)         #selfing coefficient 
p(1).General.seedDelay = params_in(7,1)            #seed delay
p(1).General.femaleFrac = params_in(8,1)           #fraction female
p(1).Init.Upper.seedDensity = params_in(9,1)       #initial upper bank density
p(1).Init.Lower.seedDensity = params_in(10,1)       #initial lower bank density

p(1).Init.Upper.ResAlleleFreq = params_in( 11:2:17,1 )   #initial upper resistant allele frequency
p(1).Init.Lower.ResAlleleFreq = params_in( 12:2:18,1 )   #initial lower resistant allele frequency

p(1).General.f_RO = params_in(19,1)          #latent time
p(1).General.f_CR = params_in(20,1)     #critical response
p(1).General.f_TR = params_in(21,1)     #complete resistance time

row = 22
for i = 1:p.General.nYears
    for j = 1:p.General.nCohorts
        p(1).Germination.survivalFraction(i,j) = params_in(row, 1)
        p(1).Cultivation.survivalFraction(i,j) = params_in(row, 2)
        p(1).Hand.survivalFraction(i,j) = params_in(row, 3)
        p(1).Mature.maxPlants(i,j) = params_in(row, 4)
        p(1).Mature.cropParam(i,j) = params_in(row, 5)
        p(1).Mature.weedCompetition(i,j) = params_in(row, 6)
        p(1).SeedProd.maxYield(i,j) = params_in(row, 7)
        p(1).SeedProd.cropParam(i,j) = params_in(row, 8)
        p(1).SeedProd.weedCompetition(i,j) = params_in(row, 9)
        p(1).Mutation.muR2S_locus1(i,j) = params_in(row, 10)
        p(1).Mutation.muS2R_locus1(i,j) = params_in(row, 11)
        p(1).Mutation.muR2S_locus2(i,j) = params_in(row, 12)
        p(1).Mutation.muS2R_locus2(i,j) = params_in(row, 13)
        p(1).Mutation.muR2S_locus3(i,j) = params_in(row, 14)
        p(1).Mutation.muS2R_locus3(i,j) = params_in(row, 15)
        p(1).Mutation.muR2S_locus4(i,j) = params_in(row, 16)
        p(1).Mutation.muS2R_locus4(i,j) = params_in(row, 17)
        p(1).Predation.survivalFraction(i,j) = params_in(row, 18)
        p(1).Winter.uppersurvivalFraction(i,j) = params_in(row, 19)
        p(1).Winter.lowersurvivalFraction(i,j) = params_in(row, 20)
        row = row + 1
    


## Read and store tilling data.
path = strcat('C:\Users\', username, '\Desktop\Tilling_Data.txt')
p(1).General.TillingFreq = dlmread(path)

## Read and store herbicide data.
path = strcat('C:\Users\', username, '\Desktop\Herbicide_Data.txt')
data = dlmread(path)

data_rows = size(data, 1)
row = 1
for i = {'Herbicide1', 'Herbicide2', 'Herbicide3', 'Herbicide4'}
    if row > data_rows
        break
    
    p(1).(i{1}).targetLocus = data(row, 1)
    for k = 1:p.General.nYears
        for m = 1:p.General.nCohorts
            p(1).(i{1}).efficacySS(k, m) = data(row, 2)
            p(1).(i{1}).efficacyRS(k, m) = data(row,3)
            p(1).(i{1}).efficacyRR(k, m) = data(row,4)
            row = row + 1
        
    



