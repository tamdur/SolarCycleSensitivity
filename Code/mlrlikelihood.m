%Use a version of mlrsensitivity with function inputs replacing user inputs
%in order to generate Figure 2
%
% Ted Amdur
% 10/4/20

%Outline
% -Generate an array of all predictor choices to be used
% -Determine optimal lags for each combination using HadCRUT median
% -Create a set of coefficients and predicted temperatures for each
% ensemble member
% -Determine covariance matrix used to assess relative likelihoods
% -Perform Bayesian likelihood estimation on each set
% -Aggregate and save variables to be used in Fig. 2 of manuscript

clearvars
f = loadallforcing(0); %Load available forcing records
lagFlag = 1; %Set to 1 to calculate lags for each combination, 0 to load past calculation
yrSt = 1882; yrEnd = 2005;
%Generate an array of all predictor choices to be used
choiceVec = allcomb(1:3,[1 3], 1:2,1:4);
choiceVec = [choiceVec zeros(size(choiceVec,1),2)];

%Determine optimal lags for each combination using HadCRUT median
if lagFlag
lagV = NaN(size(choiceVec));
for ii = 1:size(choiceVec,1)
    [bT,rnT,lags,TP] = calcmlrsensitivity(yrSt,yrEnd,0,1,1,choiceVec(ii,:),f);
    b(ii)=bT(2); rn(ii,:)=rnT; lagV(ii,1:length(lags)) = lags;
    
    % -Create a set of coefficients and predicted temperatures for each ensemble member
    [bE,TPE] = calcenssensitivity(yrSt,yrEnd,0,1,...
    choiceVec(ii,:), lagV(ii,:));
    ensSet(ii).bE = bE;
    ensSet(ii).TPE = TPE;
    ensSet(ii).b = bT;
    ensSet(ii).lags = lags;
    ensSet(ii).TP = TP;
    disp(['combination ' num2str(ii) ' completed'])
end
save('../Data/code_generated/likelihoodlags_151020.mat','lagV','b','rn','ensSet')
else
    load likelihoodlags_151020.mat;
end
    
%Following is from predictor_likelihood_2_20.m
%------------------------------------------------------------------
% Determine the temperature uncertainty over time from HadCRUT4, models
%------------------------------------------------------------------
blkS = 12; %Size of blocks to use in calculation
cmipU = cmip_T_uncertainty;
subN = floor(size(ensSet(1).TPE,2)/blkS);
subCMIP = floor(length(cmipU.SD)/blkS);
uCMIP = NaN(subN,1);
cliU = NaN(subCMIP,1);
for ii = 1:subCMIP
    %Calculate the standard deviation of the block
    interval = 1+(ii-1)*blkS:ii*blkS;
    TBlk = mean(cmipU.T(:,interval),2);
    TPBlk = mean(cmipU.TP(:,interval),2);
    diff = TBlk-TPBlk;
    cliU(ii) = std(diff);
end
uCMIP(1:subCMIP) = cliU;
if length(uCMIP) > subCMIP
    uCMIP(subCMIP+1:end) = mean(uCMIP(end-3:end));
end
if isnan(uCMIP(end)) %Calculate likelihood for 2005-present period
    lastVal = uCMIP(~isnan(uCMIP));
    lastVal = lastVal(end);
    uCMIP(isnan(uCMIP)) = lastVal;
end

%------------------------------------------------------------
% Perform autoregressive joint distribution analysis
%------------------------------------------------------------
%define coefficients 
%NOTE: Following two autocorrelations are determined empirically in
%getresidualautocorr.m as the median results
a=0.5638; %autocorrelation internal variability
b=0.9;  %autocorrelation from obs.
oVari = tempvariabilitycalc(yrSt,yrEnd,blkS); %K, external variability
%Modify internal and external standard deviation to account for
%autocorrelation
iVari = 0.2.^2.*(1-a.^2);
% oVariM = (nanmean(oVari).^2)*(1-b.^2);
n=yrEnd-yrSt+1;
%define covariance matrix
for i=1:n
  for j=1:n
    sigMat(i,j)=a.^abs(i-j)*iVari + b.^abs(i-j)*oVari(round((i+j)/2));
    if i == j
        sigMat(i,j) = sigMat(i,j) + oVari(j).*(1-b.^2);
    end
  end
end

%Perform Bayesian likelihood estimation on each set

T = ... %Load temperatures
    s1a1preparehadcrutmonthly('../Data/Temperature/Time_series/HadCRUT.4.6.0.0.monthly_ns_avg.txt');
    date = T{:,1};
    yrSeries = date.Year;
    series = T{:,2};
s = series(yrSeries >= yrSt & yrSeries <= yrEnd);
SD = size(s);
sA = (squeeze(nanmean(reshape(s,[12 SD(1)/12,SD(2)]),1)))';
[yrE,ens] = loadensemble;
ens = ens(yrE >=yrSt & yrE <= yrEnd,:);
SE = size(ens);
ensA = (squeeze(nanmean(reshape(ens,[12 SE(1)/12,SE(2)]),1)));
probMEns = NaN(size(ensSet,2),100); %For HadCRUT ensemble
probM = NaN(size(ensSet,2),1); %For HadCRUT median
bM = NaN(size(ensSet,2),5);
bMEns = NaN(size(ensSet,2),100,5);
for iC = 1:size(ensSet,2)
    SE = size(ensSet(iC).TPE');
    ensP = (squeeze(nanmean(reshape(ensSet(iC).TPE',[12 SE(1)/12,SE(2)]),1)));
    TP = (squeeze(nanmean(reshape(ensSet(iC).TP,[12 SD(1)/12,SD(2)]),1)))';
    for iE = 1:100
        probMEns(iC,iE) = log(mvnpdf(ensP(:,iE)-ensA(:,iE),0,sigMat));
    end
    probM(iC) = log(mvnpdf(sA-TP,0,sigMat));
    bM(iC,:) = ensSet(iC).b;
    bMEns(iC,:,:) = ensSet(iC).bE;
end

choices = choiceVec(:,1:4);
p_LR = probM(19); p_ASH = probM(47);
bLR = ensSet(19).b; bASH = ensSet(47).b;

save('../Data/code_generated/predictor_likelihoodcmip_10_15.mat','choices','probM',...
    'probMEns','p_LR','p_ASH','bM','bMEns','bLR','bASH')

function [bE,TPE] = calcenssensitivity(yrSt,yrEnd,annual,detrendTSI,...
    choices, fLags)                     
tsiChoice = choices(1);          %1. Lean (2000) 2. WLS (2005) 3. CMIP6 SOLARIS-HEPPA 4. NRLTSI2
%detrendTSI = 1;         %1 to detrend TSI over interval, 0 otherwise

anthropogenicChoice = choices(2);%1. Miller et al. (2014) 2. Dessler Forster (2018) 3. RCP 6.0 effCO2

volcanismChoice = choices(3);    %1. Sato 550nm mean aod 2. CMIP6 550 nm mean aod

ensoChoice = choices(4);         %1. Kaplan v2 Nino3.4 2. Kaplan v2 Nino3.4 with 8-year high pass
                        %3. MEI 4. ESRL Nino3.4 5. MEI using only prior DJF average
                        %NOTE: Option 5 only available for annual-resolution MLR
amoChoice = choices(5);          %1. HadISST AMO index

ipoChoice = choices(6);          %1. HadISST IPO indes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prepare temperature record and set of predictors for chosen set
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = loadallforcing(annual); %Load available forcing records
T = ...
s1a1preparehadcrutmonthly('../Data/Temperature/Time_series/HadCRUT.4.6.0.0.monthly_ns_avg.txt');
date = T{:,1};
yrSeries = date.Year;
series = T{:,2};
%Next, find indices for selected time interval and the associated time
%series of T anomalies
global s yrsI sA
yrsI = find(f(1).yr >= yrSt & f(1).yr <= yrEnd);
s = series(yrSeries >= yrSt & yrSeries <= yrEnd);
SD = size(s);
sA = (squeeze(nanmean(reshape(s,[12 SD(1)/12,SD(2)]),1)))';

[yrE,ens] = loadensemble;
ens = ens(yrE >=yrSt & yrE <= yrEnd,:);

%Build predictor array of selected predictor datasets
%TO CHANGE MAX BOUNDS OF LAG SEARCH, ALTER 'lMax' ENTRY IN EACH
%CONDITIONAL CORRESPONDING TO RESPECTIVE TYPE OF PREDICTOR
predI = 1;
predW(:,predI) = ones(length(f(1).yr),1);
if tsiChoice > 0
    predI = predI + 1;
    predW(:,predI) = f(1).preds(:,tsiChoice);
    lMax(predI-1) = 3;
    if detrendTSI
        predW(:,predI) = predW(:,predI) - nanmean(predW(:,predI));
        predW(f(1).yr >= yrSt-lMax(predI-1) & f(1).yr <= yrEnd,predI) = ...
            detrendnan(predW(f(1).yr >= yrSt-lMax(predI-1) & f(1).yr <= yrEnd,predI));
    end
end
if anthropogenicChoice > 0
    predI = predI + 1;
    predW(:,predI) = f(2).preds(:,anthropogenicChoice);
    lMax(predI-1) = 20;
end
if volcanismChoice > 0
    predI = predI + 1;
    predW(:,predI) = f(3).preds(:,volcanismChoice);
    lMax(predI-1) = 1;
end
if ensoChoice > 0
    predI = predI + 1;
    predW(:,predI) = f(4).preds(:,ensoChoice);
    lMax(predI-1) = 1;
end
if amoChoice > 0
    predI = predI + 1;
    predW(:,predI) = f(5).preds(:,amoChoice);
    lMax(predI-1) = 10;
end
if ipoChoice > 0
    predI = predI + 1;
    predW(:,predI) = f(6).preds(:,ipoChoice);
    lMax(predI-1) = 10;
end
for ii = 2:size(predW,2)
    predW(:,ii) = predW(:,ii) - nanmean(predW(:,ii));
end
%Adjust predictor to reflect chosen lags
predL = predW;
for ii = 2:size(predL,2)
    predL(fLags(ii-1)+1:end,ii) = predW(1:end-fLags(ii-1),ii);
end
pred = predL(yrsI,:);
bE = NaN(100,size(pred,2));
TPE = NaN(100,size(pred,1));
for ii = 1:size(ens,2)
    b = regress(ens(:,ii),pred);
    bE(ii,:) = b;
    TPE(ii,:) = nansum(repmat(b',[size(pred,1) 1]).* pred,2);
end
    
end
