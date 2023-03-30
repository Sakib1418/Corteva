def [metrics,data_output] = runSimulation(varargin)
## Initialization

    #Create a program status bar 
    h = waitbar(0, 'Please wait', 'Name', 'Plot progress', ...
    'CreateCancelBtn', 'setappdata(gcbf, ''canceling'', 1)')
    setappdata(h, 'canceling', 0)


    #Run Main def with no arguments once to get the params variable
    #populated
    if ~isempty(varargin)
        Params = varargin{1}
    else
        tmp = readData()
        Params = prepareParams_2(tmp)
    
    
    #Number of simulations to run
    if isinf(Params.General.fieldSize)
        nSim = 1
    else
        nSim = Params.General.nSim
    
    
    nYears = Params.General.nYears
    nCohorts = Params.General.nCohorts
    
    #Frequency values to find the corresponding resistance time for.
    #Includes three user defined values, as well as evenly spaced
    #values from 0.05 to 0.95 fraction resistant.
    f_RO = Params.General.f_RO
    f_CR = Params.General.f_CR
    f_TR = Params.General.f_TR
    f = 0.05:0.05:0.95
    f = [f_ROf_CRf_TR0.01f'0.99]
    f_len = length(f)
    
    counter = np.zeros(f_len,1)
    
    #Variables to store results
    #variable to locally store the seed bank population
    data_store = np.zeros(nCohorts*nYears,nSim+1)
    
    #variable for sing to the excel file (includes upper and lower seed
    #banks). NOTE: I retroactively added this variable. It probably would
    #be better code to have data_store and data_output be combined into a
    #single array.
    data_output = struct('upperbank_freq',np.zeros(nCohorts*nYears,nSim+1),...
                         'lowerbank_freq',np.zeros(nYears,nSim+1))
    
    #structure array for hold population metrics
    metrics = struct('IRF',np.zeros(nSim,1),'IRF_stats',0,...
                                         'probs',np.zeros(f_len,1),...
                                         't2R',Inf(f_len,nSim),...
                                         't2R_stats',Inf(f_len,5))
    
    #Make first column contain time values
    time_step = 1/nCohorts
    data_store(:,1) = 0:time_step:nYears-time_step
    data_output.upperbank_freq(:,1) = 0:time_step:nYears-time_step
    data_output.lowerbank_freq(:,1) = 0:nYears-1
    
    ## Perform nSim simulations, storing the results 
    for n = 1:nSim
        
        #Call main def
        [ Population,~ ] = Main(Params)
       
        ## All data for direct export
        #Sum over all genotypes to find the total resistant and total seed
        #population for each cohort/year
        R = sumResistant( Population.seedBank,Params.key )
        T = squeeze( sum(Population.seedBank,2) )
        R_lower = sumResistant( Population.lowerBank,Params.key )
        T_lower = squeeze( sum(Population.lowerBank,2) )
        
        #Calculate fraction resistant for each year/cohort
        f_res = R./T
        f_res_lower = (R_lower ./ T_lower)
        
        #Zero divided by zero errors produce NaN for f_res. Replace NaN by
        #0, since NaN errors come when populations get small.
        f_res( isnan(f_res) ) = 0
        f_res_lower( isnan(f_res_lower) ) = 0
                
        #Data is in (nYears x nCohorts) arrays currently. Rearrange so that it
        #is a linear array with (1 x nYears x nCohorts) elements of increasing
        #time.
        for c = 1:nCohorts
            data_store(c:nCohorts:,n+1) = f_res(:,c)
            data_output.upperbank_freq(c:nCohorts:,n+1) = f_res(:,c)            
        
        
        data_output.lowerbank_freq(:,n+1) = f_res_lower  
        
#         plot(data_store(:,1),data_store(:,n+1),'.-')
        
        #Use the trapezoid rule to "integrate" the resistance fraction
        #(built in data integration def "trapz")
        metrics.IRF(n) = trapz( data_store(:,1),data_store(:,n+1) ) / ( nCohorts*nYears*time_step )
        
        ## Determine if user defined thresholds were exceeded and at what
        #time
        [metrics.t2R(:,n),counter] = getResistanceTime(data_store(:,1),data_store(:,n+1),f,counter)

        #Update wait bar
        waitbar(n / nSim)
    
    
    #Delete wait bar upon completion of loop
    delete(h)
    
    #Convert resistance counts to a probability
    metrics.probs = counter / nSim
        
    ## Compute statistics
    #compute min, 10th,50th,90th percentile, and max time to resistance for
    #each resistance level
    metrics.IRF_stats = mean(metrics.IRF)
    
        for k = 1:length(f)
            #find only simulations where resistance developed
            ind = find(~isinf(metrics.t2R(k,:)))
            tmp_data = metrics.t2R(k,ind)
            
            if ~isempty( tmp_data )
                #min
                metrics.t2R_stats(k,1) = min( tmp_data )

                #10th percentile
                metrics.t2R_stats(k,2) = prctile( tmp_data,10 )

                #50th percentile
                metrics.t2R_stats(k,3) = prctile( tmp_data,50 )

                #90th percentile
                metrics.t2R_stats(k,4) = prctile( tmp_data,90 )

                #max
                metrics.t2R_stats(k,5) = max( tmp_data )
            

        
    
    
    
    ## Write data to Excel file
    
    #Get reference to open Excel sheet (still open from parameter entry) 
    h = actxGetRunningServer('Excel.Application')
    myBook = h.Workbooks.Item('WeedResistance_ParameterEntryInterface.xlsm')
    
    #Write metrics to "Results" tab
    mySheet = myBook.Sheets.Item('Results')
    
    #Clear old data (Excel is already doing this, so it is commented out
    #here)
#     mySheet.Range('F7:F13').Value = ''
#     mySheet.Range('C17:D1017').Value = ''
#     mySheet.Range('H8:M1017').Value = ''
#     
#     #Number of simulations
#     mySheet.Range('F7').Value = nSim
    
    #Probability of exceeding resistances thresholds
    mySheet.Range('F8').Value = metrics.probs(1)
    mySheet.Range('F9').Value = metrics.probs(2)
    mySheet.Range('F10').Value = metrics.probs(3)
    
    #Median times to threshold
    mySheet.Range('F11').Value = metrics.t2R_stats(1,3)
    mySheet.Range('F12').Value = metrics.t2R_stats(2,3)
    mySheet.Range('F13').Value = metrics.t2R_stats(3,3)
    
    #Average integrated risk
    mySheet.Range('F14').Value = mean(metrics.IRF)
    
    #Individual simulation integrated risk
    n = 1:nSim
    mySheet.Range( strcat('I18:I',num2str(17+nSim)) ).Value = n'
    mySheet.Range( strcat('J18:J',num2str(17+nSim)) ).Value = metrics.IRF
    
    #Time-to-resistance statistics
    mySheet.Range( strcat('B18:B',num2str(14+length(f))) ).Value = f(4:)
    mySheet.Range( strcat('C18:G',num2str(14+length(f))) ).Value = metrics.t2R_stats(4:,1:5)
    
    #Write raw data (upper and lower seed bank resistant fraction) to text
    #file
    dlmwrite('UpperBank.csv',data_output.upperbank_freq,',')
    dlmwrite('LowerBank.csv',data_output.lowerbank_freq,',')
    


def [tR,counter] = getResistanceTime(time,res_frac,f,varargin)
    #time is a (Mx1) array. res_frac is an (Mx1) array. time and res_frac
    #must have the same number of rows, otherwise an error is returned. f
    #is the resistant fraction whose time-to-resistance will be found. f
    #may be a (Px1) array. Function returns the time-to-resistance, tR,
    #which is a (Px1) array. varargin is a Px1 array which may be used to
    #keep "count" of how many times resistance developed for each element
    #of f. If resistance is found for current value of f, a counter
    #variable is incremented by 1. This is then returned in counter. This
    #is useful for counting how many times resistance developed over many
    #def calls
    
    if length(time) ~= length(res_frac)
        error('time and res_frac must have equal number of rows')
        
    else
        #Counter variable
        if ~isempty(varargin)
            counter = varargin{1}
        else
            counter = np.zeros(length(f),1)
        
        
        #number of evaluation points, P
        P = length(f)
        
        #output variable, a (Px1) array
        tR = np.zeros(P,1)
        
        #Loop through f and find corresponding times.
        for j = 1:P
            tmp = find(res_frac>f(j),1,'first')

            if tmp == 1
                #The initial resistance exceeded the threshold
                tR(j) = 0
                
                #Increment counter
                counter(j) = counter(j) + 1
            elseif ~isempty(tmp)
                #If the threshold is exceeded, interpolate to find the time

                #Temporary time and resistant fraction variables for
                #readability
                t2 = time(tmp)
                t1 = time(tmp-1)
                f2 = res_frac(tmp)
                f1 = res_frac(tmp-1)

                #Linear interpolation
                tR(j) = (t2 - t1) / (f2 - f1) * (f(j) - f1) + t1
                
                #Increment counter
                counter(j) = counter(j) + 1
            else
                #The threshold was never exceeded
                tR(j) = Inf
                
            
            
        
        
    

