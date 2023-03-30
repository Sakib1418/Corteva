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
                
            
            
        
        
    
