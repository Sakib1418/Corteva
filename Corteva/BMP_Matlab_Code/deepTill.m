function [upperBank,lowerBank] = deepTill(upperBank,lowerBank,tillFrequency,currentYear)
    if tillFrequency(currentYear) == 1
        tmp = upperBank;
        upperBank = lowerBank;
        lowerBank = tmp;
    else
        return
    end
end
	