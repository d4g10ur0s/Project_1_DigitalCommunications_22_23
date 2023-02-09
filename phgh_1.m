function [x,y1 ,y2] = phgh_1(L)

    x = randn(L,1);
    a11 = 0.9;
    a12 = 0.01;
    b = 1 ;
    y1 = filter(b, [1; -a11;],x);
    y2 = filter(b, [1; -a12;],x);

end
