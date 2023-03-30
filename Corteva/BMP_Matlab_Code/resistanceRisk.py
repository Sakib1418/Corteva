def [int_risk,res_risk,res_data,RoverT] = resistanceRisk(nSim)

#Runs the Main weed resistance code nSim times, determining the resistance
#risk. The time to resistance is found for every value specified in p_eval,
#which must be numbers between 0 and 1 (time is linearly interpolated
#between years). The def returns the population for each of nSim
#simulations in res_frac, an array of resistance propability co


## Run simulation nSim times

#Life stage to consider (seed bank)
lifestage = 'seedBank'

for n = 1:nSim
    
    #Call Main resistance def.
    #(Params stays the same for all sims, so only store it once)
    if n == 1
        [Population,Params] = Main_2()
        RoverT = np.zeros(Params.General.nYears,nSim)
        res_data = np.zeros(3,nSim)
        res_risk = np.zeros(Params.General.nYears,1)
        years = 0:Params.General.nYears-1
        
        
    else
        [Population,~] = Main_2(Params)
        
    
    
    #Count the resistants in each year, and the total population in each
    #year.
    
    #resistant population
    res_pop = sumResistant(Population.(lifestage)(:,:,),Params.key)
    
    #total population
    tot_pop = sum(Population.(lifestage)(:,:,),2)
    
    #resistant fraction stored in res_frac
    RoverT(:,n) = res_pop./tot_pop
    
    #Store resistance metrics: t_latent, t_crit, and t_total
    res_data(:,n) = Population.Res
    
    #Find the resistance risk as a def of time
    ind = find( RoverT(:,n) >= Params.General.criticalResTime )
    res_risk(ind) = res_risk(ind) + 1
    
    
    


#Convert res_risk to a fraction
res_risk = res_risk / nSim

#Find integrated res_risk
int_risk = trapz(years,res_risk)/length(years)

#Plot the results
scrz = get(0,'ScreenSize')
figure('Position',[scrz(1)+100,scrz(2)+100,scrz(3)*0.6,scrz(4)*0.6])

#All simulations
ax1 = axes('Position',[0.1,0.1,0.35,0.85])
hold all
plot(years,RoverT)
axis([0,years(),0,1])

#Resistance emergence
ax1 = axes('Position',[0.6,0.1,0.35,0.85])
t = linspace(years(1),years(),1000)
plot(t,pchip(years,res_risk*100,t),'k','LineWidth',2)
axis([0,years(),0,100])



## Subdefs
def Rpop = sumResistant(pop,key)

    #Determine number of gene loci and genotypes from key variable
    [nLoci,nGeno] = size(key)
    
    #Determine number of years and cohorts from population variable
    [nYears,~,nCohorts] = size(pop)
    
    #Initialize output array to 0
    Rpop = np.zeros(nYears,nCohorts)
    
    #Sum all with resistant allele
    for G = 1:nGeno
        for L = 1:nLoci
            if key(L,G) == 2 || key(L,G) == 1
                Rpop(:,:) = Rpop(:,:) + squeeze(pop(:,G,:))
                break
            
        
    
    

