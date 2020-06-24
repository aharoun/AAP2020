function dist = calculateDistOnQuantile(qL,qH,pdf,cdf)
    dist = zeros(length(pdf),1);

    indxL = find(cdf - qL>=0,1);
    indxH = find(cdf - qH>=0,1);
    
    if isempty(indxL)
        indxL = length(pdf); 
    end

    
    if isempty(indxH)
        indxH = length(pdf); 
    end
    
    if indxL == indxH
        dist(indxL) = pdf(indxL);
    else
        dist(indxL) = cdf(indxL) - qL;

        for i = (indxL+1):(indxH-1)
            dist(i) = pdf(i);
        end
        dist(indxH) = qH - cdf(indxH-1);
    end
end





