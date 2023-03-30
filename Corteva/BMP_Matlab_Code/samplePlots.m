function [ ] = samplePlots( )
%% Initialization

P = readData();

%create figure
figure('Position',get(0,'ScreenSize'))

%choose lifestage to plot
lifestage = 'seedBank';

h = waitbar(0, 'Please wait', 'Name', 'Plot progress', ...
    'CreateCancelBtn', 'setappdata(gcbf, ''canceling'', 1)');
setappdata(h, 'canceling', 0)

%run simulation number of times specified and plot emergence curve for each
subplot(2,2,1)
hold all
nSim = P.General.nSim;      %number of simulations to run
nRes = 0;       %number of times resistance develops after nYears


%% Run simulation specified number of times 

%Resistant fraction evaulation points
fEval = 0.05:0.05:0.95;
yearEval = zeros(nSim,length(fEval));

%Default threshold for when population is considered resistant
toler = P.General.criticalResTime;

%Create array to store resistance metric data for all simulations
resData = zeros(nSim, 3);

%Run simulation
for n = 1:nSim

    if getappdata(h, 'canceling')
        break
    end
    
    [Population,Params] = Main_2();
    plotPopulation(Population.(lifestage),Params.key,'fracR');
    
    %Compute resistance risk: fraction of times resistant population is >= 25% of
    %total after N years, (defined by nYears in defineParams() function)
    R = sumResistant(Population.(lifestage)(:,:,end),Params.key); %resistant population
    T = sum(Population.(lifestage)(:,:,end),2);    %total population (at end of year)
    
    if R(end)/T(end) >= 0.25
        nRes = nRes + 1;
    end

    %Create variable for storing the resistance risk using the Population
    %parameter (first time through loop only)
    if n==1
        ResistanceRisk = zeros(1,Params.General.nYears);
    end

    %count the resistants in each year and divide by total population
    resPop = sumResistant(Population.(lifestage)(:,:,end),Params.key);
    totPop = sum(Population.(lifestage)(:,:,end),2);
    resFrac = resPop./totPop;

    %find the years when the resistant fraction exceeds the tolerance
    ind = find(resFrac>=toler);
    ResistanceRisk(ind) = ResistanceRisk(ind) + 1;

    %find fractional year values corresponding to the
    %resistant_fraction array specified by the user
    yearEval(n,:) = sigmoidMetrics(0:length(resFrac)-1,resFrac,fEval);

    %Store resistance metric data for this simulation
    resData(n, :) = Population.Res;
    
    %Update the progress bar
    waitbar(n / nSim)
end

delete(h)

%% Store resistance metric data in Excel

range = strcat('D17:F', num2str(17+nSim-1));
h = actxGetRunningServer('Excel.Application');
myBook = h.Workbooks.Item('WeedResistance_ParameterEntryInterface.xlsm');
mySheet = myBook.Sheets.Item('Results');
mySheet.Range(range).Value = resData;

%% Add title and resistance risk text to first plot
title(strcat('Fraction resistant: ',num2str(nSim), ' simulations'))
pos = Params.General.nYears / 2;
text(pos,0.2,strcat('Resistance risk = ',num2str(nRes/nSim)));

%% Plots with +- one standard deviation window

%convert ResistanceRisk to a probability
ResistanceRisk = ResistanceRisk/nSim;

%eliminate curves where didn't or partially developed
k = 1;
while k <= size(yearEval,1)
    if sum(isnan(yearEval(k,:)))~=0
        yearEval(k,:) = [];
    else
        k = k + 1;
    end
end

%Take average, standard deviation, and store resistant_fraction in array
%meanCurve, for output
meanCurve(:,1) = mean(yearEval, 1);
meanCurve(:,2) = std(yearEval, 0, 1);
meanCurve(:,3) = fEval;

%Fit a sigmoid to the output (using subfunction sigmoidFun, below)
p = lsqcurvefit(@sigmoidFun,[0.5,5],meanCurve(:,1),meanCurve(:,3));

%Use fit parameters to create a new population variable pop, which may
%be smoother than the original simulation data.
t = linspace(0,Params.General.nYears-1);
pop = sigmoidFun(p,t);
t = t';
pop = pop';

% ax1 = axes('Position',[0.1,0.1,0.38,0.85]);
subplot(2,2,2)
hold all
% set(0,'DefaultAxesColorOrder',[0 0 1;1 0 0;0 0.4 0;1 0.6 0])

%lower window
p = lsqcurvefit(@sigmoidFun,[0.5,5],meanCurve(:,1)-2*meanCurve(:,2),meanCurve(:,3));
wind = sigmoidFun(p,t);
area(t,wind,'FaceColor',[0.7,0.7,0.7],'LineStyle','none')
plot(t,wind,'k','LineWidth',1,'LineStyle',':')

%upper window
p = lsqcurvefit(@sigmoidFun,[0.5,5],meanCurve(:,1)+2*meanCurve(:,2),meanCurve(:,3));
wind = sigmoidFun(p,t);
area(t,wind,'FaceColor',[1,1,1],'LineStyle','none')
plot(t,wind,'k','LineWidth',1,'LineStyle',':')

plot(t,pop,'k-','LineWidth',3)
%plot(t,-2*wind,'LineWidth',3,'LineStyle',':')
xlabel('Year','Interpreter','latex')
ylabel('Fraction resistant','Interpreter','latex')
ax = gca;
ax.FontSize = 16;

%% Resistance risk versus time plot

%Make data look smooth using Piecewise Cubic Hermite Interpolating
%Polynomial (PCHIP)

%Years in original simulation
n = 0:length(ResistanceRisk)-1;

%Points to do the PCHIP through
t = linspace(0,max(n));

%Call pchip
pop = pchip(n,ResistanceRisk,t);

%Plot results
t = t';
pop = pop';

% ax2 = axes('Position',[0.6,0.1,0.38,0.85]);
subplot(2,2,3)
hold all
plot(t,pop*100,'k-','LineWidth',3,'MarkerSize',20)
xlabel('Year','Interpreter','latex')
ylabel('Resistance risk $\%$','Interpreter','latex')
ax = gca;
ax.FontSize = 16;

end

function y = sigmoidFun(x,xdata)
    y = 1 ./ (1 + exp(-x(1).*xdata + x(2)));
end

function Rpop = sumResistant(pop,key)

    %Determine number of gene loci and genotypes from key variable
    [nLoci,nGeno] = size(key);
    
    %Determine number of years and cohorts from population variable
    [nYears,~,nCohorts] = size(pop);
    
    %Initialize output array to 0
    Rpop = zeros(nYears,nCohorts);
    
    %Sum all with resistant allele
    for G = 1:nGeno
        for L = 1:nLoci
            if key(L,G) == 2 || key(L,G) == 1
                Rpop(:,:) = Rpop(:,:) + squeeze(pop(:,G,:));
                break
            end
        end
    end
    
end

function nEval = sigmoidMetrics(n,f,fEval)
    
    %Number of evaluation points
    len = length(fEval);

    %Array to store output
    nEval = zeros(1,len);
    
    %Loop through evaluation points and populate nEval
    for k = 1:len
        
        %Find the first value >= current evaluation point
        ind = find(f>=fEval(k),1,'first');
        
        %If no value found or if the population initially exceeds
        %tolerance, store NaN, otherwise, find the time corresponding to
        %the current evaulation point (using linear interpolation)
        if isempty(ind)
            nEval(k) = NaN;
        elseif ind == 1
            nEval(k) = NaN;
        else
            nEval(k) = (n(ind)-n(ind-1))/(f(ind)-f(ind-1))*(fEval(k)-f(ind-1)) + n(ind-1);
        end
        
    end
end

