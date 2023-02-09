[x, y1, y2] = phgh_1(10000);

j=1;
i=1;
sqnr = [];
while i<2^10
    i = i * 2;
    disp(i);
    [xq, sqnr(j,1)] = aDeltaModulation(y2 , i);
    sqnr(j,2) = i;
    j = j + 1;
end

hold on
scatter(sqnr(:,2) , sqnr(:,1));
hold off
