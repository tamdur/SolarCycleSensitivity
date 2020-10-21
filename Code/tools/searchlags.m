function [fLags,rFine] = searchlags(lMax,predW,yrsI,s,annual)
%SEARCHLAGS Summary of this function goes here
%   Detailed explanation goes here
    if ~annual
        lMax = lMax *12; 
    end
    lInt = ceil(lMax./6); %Coarseness of original search (7 intervals for each predictor)    
    %Make array of coarse lag permutations to be tried
    lP = {};
    for ii = 1:length(lMax)
        lP = [lP; {0:lInt(ii):lMax(ii)}];
    end
    lA = allcomb(lP{:});
    [oLCoarse, rCoarse] = findlags(predW,lA,0,yrsI,s);
    %Use lags from coarse search to perform a finer-resolution search
    oLStart = oLCoarse' - ceil(lInt ./2); oLStart(oLStart < 0) = 0;
    lPF = {};
    for ii = 1:length(lMax)
        lPF = [lPF; {oLStart(ii):ceil(lMax(ii)./36):oLCoarse(ii)+ceil(lInt ./ 2)}];
    end
    lAF = allcomb(lPF{:});
    [fLags,rFine] = findlags(predW,lAF,0,yrsI,s);
end

function [oL,r2] = findlags(predW,lA,r2T,yrsI,s)
%predW original predictor array, lA array of NXM lag combinations where N
%is the number of combinations and M is the number of predictors, and r2T
%is the r2 threshold necessary to change to a new set of lags
if nargin < 3
    r2T = 0;
end
for ii = 1:size(lA,1)
    pT = predW;
    for iP = 1:size(lA,2)
        pT(lA(ii,iP)+1:end,iP+1) = pT(1:end-lA(ii,iP),iP+1);
    end
    [~,~,~,~,stats] = regress(s,pT(yrsI,:));
    if stats(1) >= r2T
        oL = lA(ii,:);
        r2T = stats(1);
    end
end
r2 = r2T;
end

