function runSimulation(varargin)
% Function simulates the weed life-cycle for specified number of
% simulations. Function reads parameter input Excel sheet in the same
% directory. The function writes simulation results to the same
% spreadsheet, as well as writing csv files containing population data. The
% user may pass the Params structure from a previous simulation, if
% desired.

%% Initialization & preliminary calculations

    %Create a program status bar to alert user how far along the simulation
    %is
    h = waitbar(0, 'Please wait', 'Name', 'Plot progress', ...
    'CreateCancelBtn', 'setappdata(gcbf, ''canceling'', 1)');
    setappdata(h, 'canceling', 0)

    %Run Main function with no arguments once to get the Params variable
    %populated, unless user has included Params as an argument of the
    %runSimulation function (in varargin).
    if ~isempty(varargin)
        Params = varargin{1};
    else
        tmp = readData();
        Params = prepareParams_2(tmp);
    end
    
    %Store the number of simulations desired in local variable. If the
    %fieldsize is infinite, simulation is deterministic, so nSim is set to
    %1.
    if isinf(Params.General.fieldSize)
        nSim = 1;
    else
        nSim = Params.General.nSim;
    end
    
    %Store the number of years and cohorts in local variables
    nYears = Params.General.nYears;
    nCohorts = Params.General.nCohorts;
    
    %Frequency values to find the corresponding resistance time for.
    %Includes three user defined values.
    f_RO = Params.General.f_RO;     %resistance onset frequency
    f_CR = Params.General.f_CR;     %critical resistance frequency
    f_TR = Params.General.f_TR;     %total resistance frequency
    f = [f_RO;f_CR;f_TR];           %array containing 3 values
    f_len = length(f);              %length of array
    
    %Counter variable for keeping track of the number of times resistance
    %thresholds are exceeded. 
    counter = zeros(f_len,1);
    
    %Variables to store results
    
    %Variable to locally store the seed bank population
    data_store = zeros(nCohorts*nYears,nSim+1);
    
    %Variables for sending to the csv files (includes upper and lower seed
    %banks). Both resistance frequency and total density are stored and
    %written to files. NOTE: I retroactively added this variable. It
    %probably would be better code to have data_store and freq_output be
    %combined into a single array.
    freq_output = struct('upperbank_freq',zeros(nCohorts*nYears,nSim+1),...
                         'lowerbank_freq',zeros(nYears,nSim+1));
    dens_output = struct('upperbank_freq',zeros(nCohorts*nYears,nSim+1),...
                         'lowerbank_freq',zeros(nYears,nSim+1));
    
    %Make first columns contain time values
    time_step = 1/nCohorts;
    data_store(:,1) = 0:time_step:nYears-time_step;
    freq_output.upperbank_freq(:,1) = 0:time_step:nYears-time_step;
    freq_output.lowerbank_freq(:,1) = 0:nYears-1;
    dens_output.upperbank_freq(:,1) = 0:time_step:nYears-time_step;
    dens_output.lowerbank_freq(:,1) = 0:nYears-1;
    
    %Structure array for holding population metrics
    metrics = struct('IRF',zeros(nSim,1),'IRF_stats',0,...
                                         'probs',zeros(f_len,1),...
                                         't2R',Inf(f_len,nSim),...
                                         't2R_stats',Inf(f_len,5),...
                                         'yearly_freq',{cell(nYears*nCohorts,6)} );
    

    
%% Perform nSim simulations, storing the results 
    for n = 1:nSim
        
        %Call main function and store Population structure array
        [ Population,~ ] = Main(Params);
       
    %% Compute metrics
        %Sum over all genotypes to find the total resistant and total seed
        %population for each cohort/year in the upper and lower banks
        R = sumResistant( Population.seedBank,Params.key );
        T = squeeze( sum(Population.seedBank,2) );
        R_lower = sumResistant( Population.lowerBank,Params.key );
        T_lower = squeeze( sum(Population.lowerBank,2) );
        
        %Calculate fraction resistant for each year/cohort
        f_res = R./T;
        f_res_lower = (R_lower ./ T_lower);
        
        %0/0 division errors produce NaN for f_res. Replace NaN by
        %0, since NaN errors come when populations approach 0.
        f_res( isnan(f_res) ) = 0;
        f_res_lower( isnan(f_res_lower) ) = 0;
                
        %Data is in (nYears x nCohorts) arrays currently. Rearrange so that it
        %is a linear array with (1 x nYears x nCohorts) elements of increasing
        %time.
        for c = 1:nCohorts
            data_store(c:nCohorts:end,n+1) = f_res(:,c);
            freq_output.upperbank_freq(c:nCohorts:end,n+1) = f_res(:,c);    
            dens_output.upperbank_freq(c:nCohorts:end,n+1) = T(:,c);  
        end
        
        %Store lower bank data
        freq_output.lowerbank_freq(:,n+1) = f_res_lower;  
        dens_output.lowerbank_freq(:,n+1) = T_lower;  
                
        %Use the trapezoid rule to "integrate" the resistance fraction
        %(built in data integration function "trapz")
        metrics.IRF(n) = trapz( data_store(:,1),data_store(:,n+1) ) / ( nCohorts*(nYears-1)*time_step );
        
        %% Determine if user defined thresholds were exceeded and at what
        %time
        [metrics.t2R(:,n),counter] = getResistanceTime(data_store(:,1),data_store(:,n+1),f,counter);

        %Update wait bar
        waitbar(n / nSim)
    end

%% Post-processing of simulation results
    %Delete wait bar upon completion of loop
    delete(h)
    
    %Convert resistance counts to a probability
    metrics.probs = counter / nSim;
        
    %% Compute statistics
    
    %Find the mean of the integrated resistance fraction from each
    %simulation
    metrics.IRF_stats = mean(metrics.IRF);
    
    %Compute min, 10th,50th,90th percentile, and max times required to
    %achieve the user specified resistance thresholds    
        for k = 1:length(f)
            %find only simulations where resistance developed
            ind = find(~isinf(metrics.t2R(k,:)));
            tmp_data = metrics.t2R(k,ind);
            
            if ~isempty( tmp_data )
                %min
                metrics.t2R_stats(k,1) = min( tmp_data );

                %10th percentile
                metrics.t2R_stats(k,2) = prctile( tmp_data,10 );

                %50th percentile
                metrics.t2R_stats(k,3) = prctile( tmp_data,50 );

                %90th percentile
                metrics.t2R_stats(k,4) = prctile( tmp_data,90 );

                %max
                metrics.t2R_stats(k,5) = max( tmp_data );
            end

        end
        
    %Compute the min, 10th, 50th, 90th percentile, and max values of
    %resistance for each year/cohort. If resistance values are all zeros,
    %this cohort/year is ignored (and the value is left blank but not 0).
	metrics.yearly_freq(:,1) = num2cell( data_store(:,1)+1 );
        for y = 1:nCohorts*nYears
            %
            ind = find(data_store(y,2:end)~=0);
            if ~isempty(ind)
                metrics.yearly_freq{y,2} = min(data_store(y,ind+1),[],2);
                metrics.yearly_freq{y,3} = prctile(data_store(y,ind+1),10,2);
                metrics.yearly_freq{y,4} = prctile(data_store(y,ind+1),50,2);
                metrics.yearly_freq{y,5} = prctile(data_store(y,ind+1),90,2);
                metrics.yearly_freq{y,6} = max(data_store(y,ind+1),[],2);
            end
        end
    
    %% Write data to Excel file
    
    %Get reference to open Excel sheet (still open from parameter entry) 
    h = actxGetRunningServer('Excel.Application');
    myBook = h.Workbooks.Item('WeedResistanceGUI.xlsm');
    
    %Write metrics to "Results" tab
    mySheet = myBook.Sheets.Item('Results');
    
    %Clear old data (Excel is already doing this, so it is commented out
    %here)
%     mySheet.Range('F7:F13').Value = '';
%     mySheet.Range('C17:D1017').Value = '';
%     mySheet.Range('H8:M1017').Value = '';
%     
%     %Number of simulations
%     mySheet.Range('F7').Value = nSim;
    
    %Probability of exceeding resistances thresholds
    mySheet.Range('F8').Value = metrics.probs(1);
    mySheet.Range('F9').Value = metrics.probs(2);
    mySheet.Range('F10').Value = metrics.probs(3);
    
    %Median times to threshold (convert to string in case value = Inf)
    mySheet.Range('F11').Value = num2str( metrics.t2R_stats(1,3) );
    mySheet.Range('F12').Value = num2str( metrics.t2R_stats(2,3) );
    mySheet.Range('F13').Value = num2str( metrics.t2R_stats(3,3) );
    
    %Average integrated risk
    mySheet.Range('F14').Value = mean(metrics.IRF);
    
    %Individual simulation integrated risk
    n = 1:nSim;
    mySheet.Range( strcat('I18:I',num2str(17+nSim)) ).Value = n';
    mySheet.Range( strcat('J18:J',num2str(17+nSim)) ).Value = metrics.IRF;
    
    %Resistance versus time statistics
    [rows,~] = size(metrics.yearly_freq);
    mySheet.Range( strcat('B18:B',num2str(17+rows)) ).Value = metrics.yearly_freq(:,1);
    mySheet.Range( strcat('C18:G',num2str(17+rows)) ).Value = metrics.yearly_freq(:,2:end);
    
    %Write raw data (upper and lower seed bank resistant fraction) to text
    %file
    dlmwrite('UpperBankFreq.csv',freq_output.upperbank_freq,',');
    dlmwrite('LowerBankFreq.csv',freq_output.lowerbank_freq,',');
    dlmwrite('UpperBankDens.csv',dens_output.upperbank_freq,',');
    dlmwrite('LowerBankDens.csv',dens_output.lowerbank_freq,',');
    
end

function [tR,counter] = getResistanceTime(time,res_frac,f,varargin)
    %time is a (Mx1) array. res_frac is an (Mx1) array. time and res_frac
    %must have the same number of rows, otherwise an error is returned. f
    %is the resistant fraction whose time-to-resistance will be found. f
    %may be a (Px1) array. Function returns the time-to-resistance, tR,
    %which is a (Px1) array. varargin is a Px1 array which may be used to
    %keep "count" of how many times resistance developed for each element
    %of f. If resistance is found for current value of f, a counter
    %variable is incremented by 1. This is then returned in counter. This
    %is useful for counting how many times resistance developed over many
    %function calls
    
    if length(time) ~= length(res_frac)
        error('time and res_frac must have equal number of rows');
        
    else
        %Counter variable
        if ~isempty(varargin)
            counter = varargin{1};
        else
            counter = zeros(length(f),1);
        end
        
        %number of evaluation points, P
        P = length(f);
        
        %output variable, a (Px1) array
        tR = zeros(P,1);
        
        %Loop through f and find corresponding times.
        for j = 1:P
            tmp = find(res_frac>f(j),1,'first');

            if tmp == 1
                %The initial resistance exceeded the threshold
                tR(j) = 0;
                
                %Increment counter
                counter(j) = counter(j) + 1;
            elseif ~isempty(tmp)
                %If the threshold is exceeded, interpolate to find the time

                %Temporary time and resistant fraction variables for
                %readability
                t2 = time(tmp);
                t1 = time(tmp-1);
                f2 = res_frac(tmp);
                f1 = res_frac(tmp-1);

                %Linear interpolation
                tR(j) = (t2 - t1) / (f2 - f1) * (f(j) - f1) + t1;
                
                %Increment counter
                counter(j) = counter(j) + 1;
            else
                %The threshold was never exceeded
                tR(j) = Inf;
                
            end
            
        end
        
    end
end
