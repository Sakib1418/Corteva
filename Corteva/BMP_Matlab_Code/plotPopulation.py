def plt_han = plotPopulation(LifeStage,key,varargin)
#Accepts a lifestage of the population structure and plots the population
#versus time. Required inputs are LifeStage (e.g. Population.germination)
#and the genotype key (e.g. Params.key). Optional argument is the plot
#type: 
#'total' = total population density 
#'totalR' = total population #possessing resistant allele (including
#heterozygotes) 
#'totalRR' = total #population possessing 2 resistant alleles 
#'fracR' = fraction of population containing resistant allele (including
#heterozygotes)
#'fracRR' = fraction of population containing 2 resistant alleles

if isempty(varargin)
    keyword = 'total'
else
    keyword = varargin{1}


[nYears,nGeno,nCohorts] = size(LifeStage)
N = nYears*nCohorts

newPop = np.zeros(N,nGeno)

k = 1
for y = 1:nYears
    for c = 1:nCohorts
        newPop(k,:) = LifeStage(y,:,c)
        k = k + 1
    


switch keyword
    case 'total'
        pltPop = sum(newPop,2)
        time = 0:N-1
        time = time / nCohorts
        ylab = 'Population, ($pl/m^2$)'
        ymax = Inf
    case 'totalR'
        pltPop = sumResistant(newPop,key)
        time = 0:N-1
        time = time / nCohorts
        ylab = 'Population, ($pl/m^2$)'
        ymax = Inf
    case 'totalRR'
        pltPop = sumResistant(newPop,key,'homozygous')
        time = 0:N-1
        time = time / nCohorts
        ylab = 'Population, ($pl/m^2$)'
        ymax = Inf
    case 'fracR'
        tmp = sumResistant(newPop,key)
        pltPop = tmp./sum(newPop,2)
        time = 0:N-1
        time = time / nCohorts
        ylab = 'Fraction resistant'
        ymax = 1
    case 'yearlyfracR'
        newPop = newPop(nCohorts:nCohorts:,:)
        tmp = sumResistant(newPop,key)
        pltPop = tmp./sum(newPop,2)
        time = 0:nYears-1
        ylab = 'Fraction resistant'
        ymax = 1
        


hold all
plt_han = plot(time,pltPop,'.-')
xlabel('Year','FontSize',20,'Interpreter','latex')
ylabel(ylab,'FontSize',20,'Interpreter','latex')
axis([0 max(time) 0 ymax])
set(gca,'FontSize',15)
hold off