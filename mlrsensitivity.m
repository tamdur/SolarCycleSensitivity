% Script to generate MLR sensitivity estimate
%
% Ted Amdur
% 9/30/20

clearvars
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USER ENTRY PORTION: Enter conditions for MLR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yrSt = 1882;            %Start of interval, inclusive. Must be >=1870
yrEnd = 2019;           %End of interval, inclusive. Must be <= current year
annual = 1;             %Set to 1 for annual resolution, 0 for monthly
lags = 1;               %1 to dynamically fit lags, 0 to use lags entered on next line
fLags = [1;120;6;4];      %Lags (in mths or yrs) to be used for each predictor
%LR08: [1;120;6;4];     %in order [TSI;ANTH;VOL;ENSO;AMO;IPO]
%Misios16: [2;0;0;0]    %if 'lags' parameter set to 0. Add or subtract
                        %entries to match number of desired predictors.
                        
%For all following choices, enter integer corresponding to the desired
%predictor dataset. If predictor is not to be used, ENTER 0. 
TChoice = 1;            %1 HadCRUT4 (1850 to present) 2 NCEP (1948-present) 
                        %3 Berkeley Earth (1750-2015) 4 GISS(1880-present)
                        %5 HadCRUT3

tsiChoice = 2;          %1 Lean (2000) 2 WLS (2005) 3 CMIP6 SOLARIS-HEPPA 4. NRLTSI2
detrendTSI = 0;         %1 to detrend TSI over interval, 0 otherwise

anthropogenicChoice = 1;%1 Miller et al. (2014) 2 Dessler Forster (2018) 3 RCP 6.0 effCO2

volcanismChoice = 1;    %1 Sato 550nm mean aod 2. CMIP6 550 nm mean aod

ensoChoice = 2;         %1 Kaplan v2 Nino3.4    2 Kaplan v2 Nino3.4 with 8-year high pass
                        %3 MEI      4 ESRL Nino3.4      5 MEI using only prior DJF average
                        %NOTE: Option 5 only available for annual-resolution MLR
                        
amoChoice = 0;          %1 HadISST AMO index    2 ERSSTv5 AMO index

ipoChoice = 0;          %1 HadISST IPO indes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ci = 0.95; iter = 10000; %Settings for phase randomization CI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prepare temperature record and set of predictors for chosen set
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = loadallforcing(annual); %Load available forcing records
if TChoice == 1
    T = ...
    s1a1preparehadcrutmonthly('../Data/Temperature/Time_series/HadCRUT.4.6.0.0.monthly_ns_avg.txt');
    date = T{:,1};
    yrSeries = date.Year;
    series = T{:,2};
elseif TChoice == 2
    [xNCEP, yNCEP, yrNCEP,mthNCEP, varNCEP] = s1apreparencep25;
    [yrSeries,series] = takespatialavgmonthly(yrNCEP,mthNCEP,yNCEP,varNCEP, 2);
elseif TChoice == 3
    [yrSeries,series] = s1a3prepareberkeley;
elseif TChoice == 4
    [xGiss, yGiss, yrGiss, mthGiss, varGiss] = s1a4preparegiss;
    [yrSeries,series] = takespatialavgmonthly(yrGiss,mthGiss,yGiss,varGiss, 2);
elseif TChoice == 5
    s1preparehadcrut3
    yrSeries = yrH; series = TH;
else
    error('Must enter valid choice of temperature time series')
end
global s yrsI
%Next, find indices for selected time interval and the associated time
%series of T anomalies
yrsI = find(f(1).yr >= yrSt & f(1).yr <= yrEnd);
s = series(yrSeries >= yrSt & yrSeries <= yrEnd);
if annual
    SD = size(s);
    s = (squeeze(nanmean(reshape(s,[12 SD(1)/12,SD(2)]),1)))';
end
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
    if ensoChoice == 5 && annual == 0
        error('ensoChoice cannot be set to 5 for monthly-resolution MLR')
    end
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
lMax = lMax';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate optimal lags, if requested
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if lags == 1
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
    [oLCoarse, rCoarse] = findlags(predW,lA,0);
    %Use lags from coarse search to perform a finer-resolution search
    oLStart = oLCoarse' - ceil(lInt ./2); oLStart(oLStart < 0) = 0;
    lPF = {};
    for ii = 1:length(lMax)
        lPF = [lPF; {oLStart(ii):ceil(lMax(ii)./36):oLCoarse(ii)+ceil(lInt ./ 2)}];
    end
    lAF = allcomb(lPF{:});
    [fLags,rFine] = findlags(predW,lAF,rCoarse);
end
%Adjust predictor to reflect chosen lags
predL = predW;
for ii = 2:size(predL,2)
    predL(fLags(ii-1)+1:end,ii) = predW(1:end-fLags(ii-1),ii);
end
pred = predL(yrsI,:);
range = phaserandommlr(s, pred, ci,iter); %Get confidence interval
[b,~,~,~,stats] = regress(s,pred);
disp(['Global T sensitivity: ' num2str(b(2),'%.3f')  ' K/(W m^2})' ...
    char(177) num2str((range(2)-range(1))/2,'%.3f') ' K/(W m^2})'])
function [oL,r2] = findlags(predW,lA,r2T)
%predW original predictor array, lA array of NXM lag combinations where N
%is the number of combinations and M is the number of predictors, and r2T
%is the r2 threshold necessary to change to a new set of lags
global yrsI s;
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

% Code to save LR08 data for Fig. 1
% bLR = b; fLagsLR = fLags; sLR = s; predLR = pred; yrFLR = f(1).yrF(yrsI);
% save('LR08to2019.mat','bLR','fLagsLR','sLR','predLR','yrFLR')
