function splitstr = strsplit(str,delimiters)
    Nd = length(delimiters);
    places = 0;
    for d = 1:Nd
        pos = find(str == delimiters{d});
        places = [places, pos]; %#ok<AGROW>
    end
    places = sort(places);

    if places(end)<length(str)
        splitstr = cell(1,length(places));
        for p = 1:length(places)-1
            splitstr{p} = str(places(p)+1:places(p+1)-1);
        end
        splitstr{p+1} = str(places(end)+1:end);
    else
        splitstr = cell(1,length(places)+1);
        for p = 1:length(places)-1
            splitstr{p} = str(places(p)+1:places(p+1)-1);
        end
    end
end