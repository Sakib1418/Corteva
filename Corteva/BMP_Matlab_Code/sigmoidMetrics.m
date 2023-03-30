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