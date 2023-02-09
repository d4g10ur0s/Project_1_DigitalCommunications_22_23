
function [xq ,sqnr] = aDeltaModulation(sima,M)
    sqnr = 0 ;
    x = interp(sima,M);
    k=1.5;
    enc = zeros(length(x),1);
    dstep = 1;%kalo einai ke to mean gia arxh
    xq = zeros(length(x), 1);%kvantismeno shma
    
    for i=1:length(x)-1
        if x(i) - xq(i) > 0
            %ean h diafora einai megaluterh tou 0 , +1 * dstep
            xq(i+1) = dstep;
        elseif x(i) - xq(i) < 0
            %ean h diafora einai mikroterh tou 0 , -1 * dstep
            xq(i+1) = -dstep;
        else
            %ananewsh tou dstep
            if xq(i+1) == xq(i)
                dstep = dstep * k;
            else
                dstep = dstep/k;
            end

            if x(i) > 0 
                xq(i+1) = dstep;
            else
                xq(i+1) = -dstep;
            end

        end
    end
    sqnr = mean(x.^2) / mean(xq.^2);
end