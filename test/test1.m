nSims = 10000;
nSteps = 10;
xVals = zeros(nSims,1);
yVals = zeros(nSims,1);
for i = 1:nSims
    x = rand * 0.4 - 0.2;
    y = rand;
    z = rand * 0.2 - 0.1;
    for j = 1:nSteps
        gVar = 0.1 * randn;
        x = x + 0.01 * sin(y) - 0.01 * cos(y) + gVar;
        y = y + 0.1 * z + rand * 0.2 - 0.1;
        z = rand * 0.2 - 0.1;
    end
    xVals(i,1) = x;
    yVals(i,1) = y;
end
scatter(xVals,yVals,'.');
disp(mean(xVals));
disp(mean(yVals.^2));
disp(mean(yVals));
