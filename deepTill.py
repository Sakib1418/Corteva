def deepTill(upperbank,lowerbank,tillFrequency,currentYear):
    if tillFrequency[currentYear] ==1:
        tmp = upperbank
        upperbank = lowerbank
        lowerbank = tmp 
    else:
        pass
    return upperbank, lowerbank