
function [SER ,BER , out] = mpam(bslen,bitn,SNR)
    %0. Xronikes Parametroi, kanonikopoihmenes gia thn prosomoiwsh
    A =1 ;
    
    Rsymbol = 250 * 10e+3;
    Tsymbol = 40;
    Tc = 4;
    fc=Tc/4;
    Tsample = 1;

    %1. binary sequence
    bs = binary_sequence(bslen);
    %2. create gray code
    code = grayCode(bitn);
    %3. map symbols
    [frequency , mapping] = mMapper(bs, code);
    %4-5. modulation
    disp("modulation");
    %modulation
    gt = sqrt(2 / Tsymbol);%rect pulses
    sm = zeros(strlength(mapping) * 40, 1);
    t = (0:Tsample:(Tsymbol - Tsample))';
   
    j=1;
    while j < strlength(mapping)
        e = str2num(extract(mapping,j));
        for k=1:length(t)
            sm((j - 1) * 40 + k) = (2*e - 1 - log2(bitn)) * gt * cos(2 * pi * fc * t(k));    
        end
        j=j+1;
    end

    %6. awng
    Eb = 1/log2(bitn);
    var = Eb/2 * 10.^(-SNR/10);
    noise = sqrt(var) .* randn(length(sm), 1);
    
    r = sm+noise;
    
    %7. demodulation
    disp("demodulation")
    for i=1:length(r)
        r(i)=r(i)*gt*cos(2*pi*fc*(i-1)*Tsample)*Tsample;
    end

    demodulated=zeros(size(sm));
    
    for i=1:strlength(mapping)
        indx=((i-1)*40)+1;
        demodulated(i)=sum(r(indx:(indx+39)));
    end
    
    %8. decision
    disp("decision")
    decision=zeros(strlength(mapping),1);
    for i=1:strlength(mapping)
        dist=abs( (2*(0:1:log2(bitn)) - 1 - log2(bitn))-demodulated(i));
        [~,ind]=min(dist);
        decision(i)=ind; 
    end

    %9. demapping
    disp("demapping");
    outp = mDeMapper(decision,code);
    %disp(outp);
   
    % Calculate SER
    errorsym=0;
    totalsymb=0;
    
    for i=1:strlength(mapping)
        if decision(i)-str2num(extract(mapping,i))~=0
            errorsym=errorsym+1;
        end
        totalsymb=totalsymb+1;
    end
    %This is the end
    mend=length(bs);
    k=abs(sum(bs(1:mend-1)-outp) );
    
    BER=length(k)/length(bs);
    out=outp;
    SER=errorsym/totalsymb;

end

function [frequency, mapping] = mMapper(bsequence , code)
    codelen = log2(length(code));
    frequency = zeros(length(code), 1);
    mapping = "";
    %to iteration 8a ksekinhsei gia i=1 alla 8a auksanetai me mhkos kwdika
    i=1;
    while i<length(bsequence)
        j=1;
        while j<length(code)+1
            if strcmp(toString(bsequence(i:i+codelen-1)), code{j})    
                mapping = mapping + string(j);
                frequency(j,1) = frequency(j,1) + 1;
                break;
            end
            j = j+1;
        end
        i = i+codelen;
    end

end

function demmaped = mDeMapper(sm ,code)
    demmaped = zeros(strlength(code{1})*length(sm),1);
    i=1;
    while i<length(sm)
        c = code{sm(i)+1};
        for j=1:strlength(code{sm(i)+1})
            demmaped(i+j-1) = str2num(extract(c,j));
        end
        i=i+strlength(code{sm(i)+1});
    end
end

function str = toString(vec)
    str = "";
    for i=1:length(vec)
        str = str + string(vec(i));
    end
end

function bs = binary_sequence(N)

    %generate random numbers then floor
    bs = floor(unifrnd(0,2,N,1));

end


function code = grayCode(M)
    
    %1. ksekinas me 1 epipedo
    code = {"0", "1"};
    %2. Gia ka8e epipedo
    for i=1:log2(M)-1
        ncode = {};
        %3. Gia ka8e string , pare to teleutaio num
        for j=1:2^i
            the_code_1 = "";
            the_code_2 = "";
            num = extract(code{j},strlength(code{j}));
            %4. An i>2 sto endiameso paei anapoda
            if i > 2
                if j == i/2
                    if num=="1"
                        the_code_2 = the_code_2 + code{j} + "1";
                        the_code_1 = the_code_1 + code{j} + "0";
                    else
                        the_code_1 = the_code_1 + code{j} + "1";
                        the_code_2 = the_code_2 + code{j} + "0";
                    end
                else
                    if num=="1"
                        the_code_1 = the_code_1 + code{j} + "1";
                        the_code_2 = the_code_2 + code{j} + "0";
                    else
                        the_code_2 = the_code_2 + code{j} + "1";
                        the_code_1 = the_code_1 + code{j} + "0";
                    end
                end
            %4. An i<2 ola kala
            else
                if num=="1"
                    the_code_1 = the_code_1 + code{j} + "1";
                    the_code_2 = the_code_2 + code{j} + "0";
                else
                    the_code_1 = the_code_1 + code{j} + "0";
                    the_code_2 = the_code_2 + code{j} + "1";
                end
            end
            ncode{2*j-1} = the_code_1;
            ncode{2*j} = the_code_2;
        end
        code = ncode;
    end
    
end