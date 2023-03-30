function P = prepareParams_2(Params)
%Takes the user defined parameter structure array created in readData and 
%prepares three variables for the simulation: UpperSeedBank, LowerSeedBank 
%and Params. UpperSeedBank and LowerSeedBank contain the initial population
%information.

%% General parameters

P.General = Params.General;
P.Init = Params.Init;
P.Cultivation = Params.Cultivation;
P.Germination = Params.Germination;
P.Hand = Params.Hand;
P.Mature = Params.Mature;
P.SeedProd = Params.SeedProd;
P.Mutation = Params.Mutation;
P.Predation = Params.Predation;
P.Winter = Params.Winter;
P.Herbicide1 = Params.Herbicide1;
P.Herbicide2 = Params.Herbicide2;
P.Herbicide3 = Params.Herbicide3;
P.Herbicide4 = Params.Herbicide4;

%Define a few of the general parameters locally for use in
%the function (not most computationally efficient, but improves
%readability)
nLoci = Params.General.nLoci;

%Generate a key for relating each genotype to the number of resistant
%alleles it posseses (0,1, or 2) at each gene locus. (Subfunction contained
%in this file).
[P.key,P.keyText] = generateKey(nLoci);


%% Initialize the upper seed bank

%Calculate the fraction of each genotype by calling the matingEquilibrium
%function (subfunction in this file).
tmp = matingEquilibrium(P.Init.Upper.ResAlleleFreq,P.General.selfingCoeff,P.key);

%Populate UpperSeedBank 
P.General.UpperSeedBank = tmp*Params.Init.Upper.seedDensity;

%% Initialize the lower seed bank

%Calculate the fraction of each genotype by calling the matingEquilibrium
%function (subfunction in this file).
tmp = matingEquilibrium(P.Init.Lower.ResAlleleFreq,P.General.selfingCoeff,P.key);

%Populate LowerSeedBank 
P.General.LowerSeedBank = tmp*Params.Init.Lower.seedDensity;

end

%% Subfunctions
function [key,key_text] = generateKey(nLoci)
%create a "key" array with rows correpsonding to a single gene locus,
%columns corresponding to unique genotypes, and each element containing
%the number of resistant alleles that a genotype has at a given gene
%locus. Example: key = [  0,0,0,1,1,1,2,2,2 ;
%   0,1,2,0,1,2,0,1,2 ] corresponds to the 9 unique genotypes when 2
%   loci are considered. First column contains 0 resistant alleles at
%   loci 1 and loci 2.
%This function also creates a cell array which contains the same
%information in more readable text format. This way, the user can
%easily correlate columns in variable "key" to a genotype. Example:
%key_text = [S1S1S2S2, S1S1R2S2, S1S1R2R2, ...etc, just run the code to
%see the output.

%number of unique genotypes
nGeno = 3^nLoci;

%initialize key array
key = zeros(nLoci,nGeno);

%populate the key array using the function odometer
tmp = zeros(nLoci,1);
for k = 1:nGeno
    key(:,k) = tmp;
    tmp = odometer(tmp,2);
end
key_text = 0;

%create the text version for easier viewing by user
key_text = cell(1,nGeno);
for k = 1:nGeno
    tmp = '';
    for j = 1:nLoci
        locus_string = num2str(j);
        if key(j,k) == 0
            tmp = strcat(tmp,'S',locus_string,'S',locus_string);
        elseif key(j,k) == 1
            tmp = strcat(tmp,'R',locus_string,'S',locus_string);
        else
            tmp = strcat(tmp,'R',locus_string,'R',locus_string);

        end

    end

    key_text{k} = tmp;

end
end
    
function odo = odometer(odo,base_num)
%% Odometer function
%This function counts using a user specified number base. A digit takes
%on integer values between 0 and base_num. This function increments the
%right-most digit, unless the right-most digit is at the maximum value,
%in which case it is set to the first element of base_num and the next
%element is incremented. For example: odo = [0,0,1] will become odo =
%[0,0,2]. A second function call on odo will produce odo = [0,1,0].
%This works like a standard car odometer, but with an arbitrary number
%system.

N = length(base_num);       %number of possible elements
Ndigits = length(odo);      %number of elements in odo

j = 1;
while j <= Ndigits
    if odo(j) < base_num
        odo(j) = odo(j)+1;
        break
    end
    odo(j) = 0;
    j = j + 1;
end
end

function populationFraction_outCross = matingEquilibrium(p,selfCoeff,key)
q = 1-p;
nLoci = size(key,1);
nGeno = size(key,2);
populationFraction_outCross = ones(1,nGeno);
%     populationFraction_selfing = ones(1,nGeno);

for j = 1:nLoci
    for k = 1:nGeno
        if key(j,k) == 2
            populationFraction_outCross(k) = populationFraction_outCross(k)*p(j)^2;
        elseif key(j,k) == 1
            populationFraction_outCross(k) = populationFraction_outCross(k)*2*p(j)*q(j);
        else
            populationFraction_outCross(k) = populationFraction_outCross(k)*q(j)^2;
        end
    end
end
end
