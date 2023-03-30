function [ resdata ] = resMetrics( resdata, currentYear, seedPop, ...
                                        key, onset, critical, complete )
                                    
%resMetrics updates weed resistance data using population data from the 
%current year of the simulation. It takes as parameters the resistance
%metrics array, the current year, the mature population for the current
%year, the key, and the user-defined threshold values (percentages at which
%weed resistance is considered to have developed, reached a critical value,
%and affected the whole population).

R = sumResistant(seedPop(currentYear, :, end),key); %resistant population
T = sum(seedPop(currentYear, :, end),2);    %total population (at end of year)

%Checks if: 1) the fraction of the populattion that is resistant has exceeded 
%the onset value and 2) the time of resistance onset has not been set. If
%both conditions are met, the time of resistance onset is set to the
%current year.
if R(end)/T(end) >= onset && resdata(1, 1) == 0
    resdata(1, 1) = currentYear;
end

%Checks if: 1) the fraction of the populattion that is resistant has exceeded 
%the critical value and 2) the time of critical resistance has not been set. 
%If both conditions are met, the time of critical resistance is set to the
%current year.
if R(end)/T(end) >= critical && resdata(1, 2) == 0
    resdata(1, 2) = currentYear;
end

%Checks if: 1) the majority of the population is resistant (according to 
%the value specified by the user) and 2) the time of complete resistance has
%not been set. If both conditions are met, the time of complete resistance 
%is set to the current year.
if R(end)/T(end) >= complete && resdata(1, 3) == 0
    resdata(1, 3) = currentYear;
end

end

