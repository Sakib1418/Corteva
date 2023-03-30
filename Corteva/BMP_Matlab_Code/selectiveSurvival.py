def newPop = selectiveSurvival(oldPop,fSurvival,targetLocus,key)
#fSurvival = 3 x 1 array containing: [fSS,fRS,fRR]
[nLoci,nGeno] = size(key)
newPop = oldPop*0

if targetLocus > nLoci
    warning('Target locus exceeds the number specified by the user in the input file, original population returned')
    newPop = oldPop
else
   for k = 1:nGeno
       switch key(targetLocus,k)
           case 0
               #homozygous susceptible
                newPop(k) = fSurvival(1)*oldPop(k)
           case 1
               #heterozygous
               newPop(k) = fSurvival(2)*oldPop(k)
           case 2
               newPop(k) = fSurvival(3)*oldPop(k)
       
       
   
   
