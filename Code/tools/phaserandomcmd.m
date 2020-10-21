function range = phaserandomcmd(Yr, Y,yI,tsi,diffwshift,varN,vol,yrSt,yrEnd, interval,iter,yrCensor,cmdMethod)
%PHASERANDOMCMD Produce lower and upper bounds of confidence interval for
%CMD using iter instances of phase randomization
%   INPUTS:
%       Yr: NX1 Year vector
%       tsi: NX1 vector of processed TSI corresponding to Yr
%       diffwshift: Whether to apply difference with shift process
%       varN: YXN array of demeaned zonal mean temperature anomalies with N
%       columns corresponding to Yr vector
%       vol: structure containing volcanic eruption data
%       yrSt: first year of interval
%       yrEnd: last year of interval
%       interval: Size of CI. Enter 0.95 for 95%
%       iter: number of iterations of phase randomization
%   OUTPUTS:
%       range: 2X1 vector of lower bound and upper bound for CI
rng(1);
sR = NaN(iter,1);
for ii = 1:iter
    s = getcmdsensitivity(Yr, Y,yI,phaserandomize(tsi),diffwshift,varN,vol,yrSt,yrEnd,yrCensor,cmdMethod);
    sR(ii) = s;
end
range(1) = quantile(sR,0.5-interval/2);
range(2) = quantile(sR,0.5+interval/2);
end

function s = getcmdsensitivity(Yr,Y,yI, tsi,diffwshift,varN,vol,yrSt,yrEnd,yrCensor,cmdMethod)
%Remove years close to zero anomaly
if diffwshift %Make two sets of years based upon shift 
    [YrN,lowSet,highSet] = makepairwise(Yr,tsi,0);%no shift
    [YrN2,lowSet2,highSet2] = makepairwise(Yr,tsi,1);%shift
    iKeep1 = ismember(Yr,YrN); iKeep2 = ismember(Yr,YrN2);
    Yr2 = Yr(iKeep2); varN2 = varN(:,iKeep2); tsi2 = tsi(iKeep2);
    Yr = Yr(iKeep1); varN = varN(:,iKeep1); tsi = tsi(iKeep1);
else
    thresh = 0.049; %Threshold of TSI anomaly beyond which a year is counted
    Yr = Yr(abs(tsi) >= thresh);
    varN = varN(:,abs(tsi) >= thresh);
    tsi = tsi(abs(tsi) >= thresh);
    lowSet = tsi <= -thresh;
    highSet = tsi >= thresh;
end
if yrCensor == 1 %Remove two years following major volcanic eruptions 
    %Find all the volcanic years in the stYr:endYr period
    volInd = find(vol.year >= yrSt & vol.year <= yrEnd);
    volYears = vol.year(volInd);
    volVals = vol.val(volInd);
    volCut = quantile(vol.val, 0.9); %percentile of cutoff for volcanic index
    volYears = volYears(volVals<=volCut);
    volYears = [volYears+1; volYears + 2]; %SET TO REFLECT CAMP/TUNG METHOD
    iRem = ~ismember(Yr,volYears);
    if diffwshift
        iRem2 = ~ismember(Yr2,volYears);
    end
elseif yrCensor == 2 %Censor years deemed volcanic by Camp and Tung (2007) 
    volYears = [1982;1983;1992;1993];
    iRem = ~ismember(Yr,volYears);
    if diffwshift
        iRem2 = ~ismember(Yr2,volYears);
    end
else
    iRem = ones(length(Yr),1);
end
%Make two groups reflecting the + and - states of solar cycle
lVec = iRem.*lowSet; hVec = iRem.*highSet;
iKeep = logical(lVec+hVec);
Yr = Yr(iKeep); tsi = tsi(iKeep); varN = varN(:,iKeep);
G = logical([lVec(iKeep) hVec(iKeep)]);
if  diffwshift
    lVec2 = iRem2.*lowSet2; hVec2 = iRem2.*highSet2;
    iKeep2 = logical(lVec2+hVec2);
    Yr2 = Yr2(iKeep2); tsi2 = tsi2(iKeep2); varN2 = varN2(:,iKeep2);
    G2 = logical([lVec2(iKeep2) hVec2(iKeep2)]);
end

%Perform CMD (traditional or LDA)
wVec = cos((pi/180).*Y(yI));
wVec(abs(Y(yI)) == 90) = 0;
wVec = wVec./mean(wVec);

if diffwshift
   cmd1 = nanmean(varN(:,G(:,2)),2) - nanmean(varN(:,G(:,1)),2);
   cmd2 = nanmean(varN2(:,G2(:,2)),2) - nanmean(varN2(:,G2(:,1)),2);
   cmd = (cmd1 + cmd2) ./ 2;
else
   cmd = nanmean(varN(:,G(:,2)),2) - nanmean(varN(:,G(:,1)),2);
end
cmd = cmd ./ nanmean(cmd.*wVec);
if cmdMethod == 1
   regSeries = [];
   for ii = 1:length(Yr)
       yrVarW = varN(:,ii).*wVec;
       regSeries = [regSeries; nansum(yrVarW.*cmd)/nansum((cmd.^2).*wVec)];
   end
   b = regress(regSeries,[ones(length(Yr),1) tsi]);
   s = b(2);
   correlation = corr(regSeries,tsi);
elseif cmdMethod == 2 %LDA
    %Define length of observations, number of groups
    n = length(Yr); %number of observations per variable
    g = 2; %number of groups
    %Perform a centering operation
    meanVar = squeeze(nanmean(varN,2));
    meanVar = repmat(meanVar,[1,n]);
    %Finish nxp array of temperatures vX where n is obs, p is all observations
    vX = (varN - meanVar)';
        wVecM = zeros(length(Y(yI)),length(Y(yI)));
    wVecM(1:length(Y(yI))+1:end) = sqrt(wVec);
    vXs = vX*wVecM; %Area-weight temperature observation vector
    %Define matrix of group means, mu^T. This matrix is [num groups X num obs]
    M = (G' * G)^(-1) * G' * vXs;
    %Create within-group covariance matrix 
    sW = (1/(n-g))*(vXs - G*M)' * (vXs - G*M);
    %Create between-group covariance matrix
    sA = (g/(n*(g-1)))*(G*M)' * (G*M);
    %Create total-covariance matrix
    sT = (1/(n-1))*((n-g)*sW + ((n*(g-1))/g)*sA);
    %Perform SVD to get first 17 eigenvectors
    [~,~,dB] = svd(sT);
    %Transform data to this basis
    if size(dB,1) > 17
        b17 = dB(:,1:17);
    else
        b17 = dB;
    end
    vXsN = vXs*b17;
    %Define NEW matrix of group means, mu^T. This matrix is [n groups X p obs]
    MN = (G' * G)^(-1) * G' * vXsN;
    %Create NEW within-group covariance matrix 
    sWN = (1/(n-g))*(vXsN - G*MN)' * (vXsN - G*MN);
    %Create NEW between-group covariance matrix
    sAN = (g/(n*(g-1)))*(G*MN)' * (G*MN);
    %Create NEW total-covariance matrix
    sTN = (1/(n-1))*((n-g)*sWN + ((n*(g-1))/g)*sAN);
    %Perform minimization operation, starting with a guess
    [~,~,VT] = svd(sTN\sAN);
    uT = VT(:,1);
    %sepMin = -1*inv(uT'*sTN*uT)*(uT'*sAN*uT); %the -1 used to make it a minimizing function
    options = optimoptions('fminunc','Display','off');
    [uSep,gammaM] = fminunc(@(uSep)-1*(uSep'*sTN*uSep)\(uSep'*sAN*uSep),uT,options);
    u72 = b17*uSep;
    u72 = u72./ sqrt(((u72)'*sT*(u72)));
    denom72 = sum((u72).^2);
    regSeries = vXs*u72./denom72;
    %Find disriminant pattern, normalize it to 1 to get c(t) magnitude
    cmd = regSeries'*vX./(regSeries'*regSeries);
    sF = nanmean(cmd'.*wVec);
    cmd = cmd ./sF; regSeries = regSeries.*sF;
    correlation = corr(regSeries,tsi);
    b = regress(regSeries,[ones(length(Yr),1) tsi]);
    s = b(2);
end


end

function [Yr,lowSet,highSet] = makepairwise(Yr,tsi,shift)
    %Yr: Years
    %tsi: tsi timeseries corresponding to years in Yr
    %shift: whether this go-around is pairwise difference with or w/o shift
    cycSt = Yr(1); cycEnd = Yr(end);
    [setSt,setEnd] = s1c4makepairwiseyears(cycSt,cycEnd,shift);
    yrsIncluded = [];
    lowYears = []; highYears = [];
    for ii = 1:length(setSt)
        iVar = find(Yr >= setSt(ii) & Yr <= setEnd(ii));
        tsiSub = tsi(iVar); tsiSub = tsiSub - mean(tsiSub);
        subYr = Yr(iVar);
        %Select years which correspond to solar max pool, solar min pool
        highInd = tsiSub >= 0.06;
        lowInd = tsiSub <= -0.06;
        %Frame shift to match the years that were passed in 
        lowYears = [lowYears;subYr(lowInd)]; 
        highYears = [highYears;subYr(highInd)];
        yrsIncluded = [yrsIncluded; lowYears; highYears];
    end
    Yr = unique(yrsIncluded);
    lowSet = ismember(Yr,lowYears);
    highSet = ismember(Yr,highYears);
end


