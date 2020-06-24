
function [managerShareModelBySize,qq] = calculateManagerShareBySize(eq)
    % Calculate share of managers by size  and quantile of size for India
    pdf = eq.fAll./sum(eq.fAll);
    cdf = cumsum(pdf);

    % calculate manager share on a regular quantile
    qq  = [0.0 .25 .5 .75 .9 .95 .99 .999];

    for i = 1:(length(qq))
        dist = calculateDistOnQuantile(qq(i),1.00000000000000,pdf,cdf);
        managerShareModelBySize(i) = (dist'*eq.manager)/(dist'*eq.emp);
    end
end
