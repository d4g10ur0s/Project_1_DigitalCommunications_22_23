function [xq , centers , D, sqnr] = lmax(x, N , minval , maxval)
    cmean = 0;
    %arxikopoihsh
    xq=zeros(length(x),1);

    kmax = 1000;%plh8os epanalhpsewn algori8mou
    er = 10e-9;%max error
    
    M = N^2;%plh8os zwnwn kvantishs
    Tk = zeros(kmax,M);%oria zwnwn kvantishs, kratw ta oria ka8e epanalhpshs
    
    qlevels = unifrnd(minval,maxval,M,1);%epilogh tyxaiwn kentrwn
    qlevels = sort(qlevels);
    
    i=1;
    while i < kmax +1
        %1. Ypologismos zwnwn kvantishs
        Tk(i,:) = compute_limits(qlevels,minval,maxval);
        %2. Ypologismos kvantismenou shmatos
        [xq, cmean] = quantize_sig(x, Tk(i,:) , qlevels);
        %3. Ypologismos newn epipedwn kvantishs
        qlevels = new_qlevels(cmean, minval , maxval,M);
        %4. Sugklish ?
        D(i) = immse(x,xq);
        sqnr(i,1) = mean(x)^2 / mean(xq)^2;
        if i > 1
            if ( D(i) - D(i-1) )^2 < er
                break;
            end
        end

        i=i+1;%aukshsh i
    end
    
    
    % Ypologismos entropias
    p = zeros(length(qlevels),1);
    for i=1:length(x)
        %iteration over qlevels
        for j=1:length(qlevels)
            if xq(i)==qlevels(j)
                p(j) = p(j) + 1;
            end
        end
    end
    
    %an einai ola isopi8ana, tote ka8e sample exei length(qlevels)
    %epiloges, ola exoun length(qlevels)^length(x)
    p = p / ( length(x) );
    entr=0;
    for i=1:length(p)
        if (p(i) == 0)
            entr = entr + 0;
        else
            entr = entr + (-p(i) * log2( p(i) ) );
        end
    end

    fprintf('Entropy : %f', entr);

    centers = qlevels;%epistrefw ta shmeia kvantismou
end

function Tk = compute_limits(qlevels,minval,maxval)
    %arxikopoihsh
    Tk = zeros(1,length(qlevels));
    %ypologismos kentrwn
    for i = 1:length(qlevels)
        if i == 1
            Tk(1,i) = (minval + qlevels(i))/2;
        end
        if i == length(qlevels)
            Tk(1,i) = (maxval + qlevels(i))/2;
        else
            Tk(1,i) = ( qlevels(i)+qlevels(i+1) )/ 2;
        end
    end
end


function [xq , cmean] = quantize_sig(x, Tk , qlevels)
    xq = zeros(length(x),1);
    %arxikopoihsh conditional mean
    cmean = zeros(length(qlevels),1);%to x na einai sthn perioxh kvantishs
    counter = zeros(length(qlevels),1);

    %gia ka8e element sto x
    for i=1:length(x)
        %gia ka8e zwnh kvantishs
        j=1;
        while j<=length(Tk(1,:))
            %an einai mikrotero apo elaxisth timh
            if j==1 && x(i)<Tk(1,1)
                xq(i,1) = qlevels(j,1);
                cmean(j,1) = cmean(j,1) + x(i,1);
                counter(j,1) = counter(j,1) + 1;
                break;
            end
            %an einai megalutero apo thn megisth timh
            if j==length(Tk(1,:)) 
                xq(i,1) = qlevels(j,1);
                cmean(j,1) = cmean(j,1) + x(i,1);
                counter(j,1) = counter(j,1) + 1;
                break;
            end
            %stis endiameses perioxes
            if j>1
                if x(i) >= Tk(1,j-1) && x(i) < Tk(1,j)
                    xq(i,1) = qlevels(j,1);
                    cmean(j) = cmean(j) + x(i);
                    counter(j) = counter(j) + 1;
                    break;
                end
            end
            j = j+1;
        end%end while
    end%end for
    
    %prepei na upologisw conditional mean
    %yparxoun ke mhdenika opote den glitwnw iteration
    for i=1:length(qlevels)
        if counter(i) == 0
            counter(i) = 1;
        end
    end
    cmean = cmean/counter;
    cmean = cmean/length(x);
end


function qlevels = new_qlevels(cmean, minval , maxval,M)
    qlevels = zeros(M, 1);
    %ta kentroeidh 
    for i = 1:M
        if i == M 
            qlevels(i)=(maxval + cmean(i))/2;
        elseif i==1
            qlevels(i)=(minval + cmean(i))/2;
        else
            qlevels(i) = (cmean(i) + cmean(i+1))/2;
        end
    end
end

