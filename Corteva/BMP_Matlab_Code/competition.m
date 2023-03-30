function newPop = competition(oldPop,A,B,C)
%Takes a population oldPop and
%applies the function newPop = A*oldPop / (1+B+C*totalPop). If oldPop is an
%1 x nGeno array, function will return a 1 x nGeno array. If oldPop is a
%total population, newPop will return a single number.

%calculate the total population
totalPop = sum(oldPop);

%calculate the new population after competition
newPop = A*oldPop ./ (1+B+C*totalPop);