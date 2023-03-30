function Rpop = sumResistant(pop,key,varargin)

if isempty(varargin)
    sumOver = 'all';
else
    sumOver = varargin{1};
end

%determine number of gene loci and genotypes from key variable
[nLoci,nGeno] = size(key);

%determine number of years and cohorts from population variable
[nYears,~,nCohorts] = size(pop);

%initialize output array to 0
Rpop = zeros(nYears,nCohorts);

switch sumOver

    case 'homozygous'
        for G = 1:nGeno
            for L = 1:nLoci
                if key(L,G) == 2 
                    Rpop(:,:) = Rpop(:,:) + squeeze(pop(:,G,:));
                    break
                end
            end
        end
        
    otherwise
        for G = 1:nGeno
            for L = 1:nLoci
                if key(L,G) == 2 || key(L,G) == 1
                    Rpop(:,:) = Rpop(:,:) + squeeze(pop(:,G,:));
                    break
                end
            end
        end
        
end