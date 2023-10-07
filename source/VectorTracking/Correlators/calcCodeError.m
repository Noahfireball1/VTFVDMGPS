function bitError = calcCodeError(pdiTime,correlatorType,codeError,freqError,chipOffset,phase)

switch upper(correlatorType)
    case 'I'
        bitError = ((1 - abs(codeError + chipOffset))*cos(pi*freqError*pdiTime + 2*pi*phase));

    case 'Q'
        bitError = ((1 - abs(codeError + chipOffset))*sin(pi*freqError*pdiTime + 2*pi*phase));

    otherwise
        bitError = 0;
end

end