function ordinal_number = ordinal_suffix(number)
digits = convertStringsToChars(string(number));
if isscalar(digits)
    digits = ['0',digits];
end

if digits(end-1) == '1'
    suffix = "th";
    ordinal_number = string(number) + suffix;
    return
end
switch digits(end)
    case '1'
        suffix = "st";
    case '2'
        suffix = "nd";
    case '3'
        suffix = "rd";
    otherwise
        suffix = "th";
end
ordinal_number = string(number) + suffix;
end