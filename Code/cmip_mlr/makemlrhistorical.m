%Run cmipreg for all cmip historical runs, returning set of all estimated 
%sensitivities.
%
% Ted Amdur
% 10/5/20, drawing from generateregcmipmodelshist_1_2020.,

%NOTE: Change cmipreg function pred matrix to use different set of forcings

%Note that recSt and recEnd has to be between 1850,2012 to fit miller
%forcings
clearvars
yrSt = 1882; yrEnd = 2005;
annual = 0; %Perform annual vs monthly mlr on CMIP runs
tsiChoice = 2;          %1 Lean (2000) 2 WLS (2005) 3 CMIP6 SOLARIS-HEPPA 4. NRLTSI2
detrendTSI = 1;         %1 to detrend TSI over interval, 0 otherwise
anthropogenicChoice = 1;%1 Miller et al. (2014) 2 Dessler Forster (2018) 3 RCP 6.0 effCO2
volcanismChoice = 1;    %1 Sato 550nm mean aod 2. CMIP6 550 nm mean aod

%hist files are stored in:
pathsTas = dir('/Users/teda/Data/Odyssey/historical45/*.nc');
pathsTos = dir('/Users/teda/Data/Odyssey/tos45/*.nc');

f = loadallforcing(annual); %Load available forcing records
%Build predictor array of selected predictor datasets
%TO CHANGE MAX BOUNDS OF LAG SEARCH, ALTER 'lMax' ENTRY IN EACH
%CONDITIONAL CORRESPONDING TO RESPECTIVE TYPE OF PREDICTOR
predW(:,1) = ones(length(f(1).yr),1);
predW(:,2) = f(1).preds(:,tsiChoice);
lMax(1) = 3;
predW(:,2) = predW(:,2) - nanmean(predW(:,2));
predW(f(1).yr >= yrSt-lMax(1) & f(1).yr <= yrEnd,2) = ...
    detrendnan(predW(f(1).yr >= yrSt-lMax(1) & f(1).yr <= yrEnd,2));
predW(:,3) = f(2).preds(:,anthropogenicChoice);
lMax(2) = 20;
predW(:,4) = f(3).preds(:,volcanismChoice);
lMax(3) = 1;
predW(:,5) = NaN(size(predW,1),1);
lMax(4) = 1;

for ii = 2:size(predW,2)
    predW(:,ii) = predW(:,ii) - nanmean(predW(:,ii));
end
yrsI = find(f(1).yr >= yrSt & f(1).yr <= yrEnd);

for iName = 1:length(pathsTas)
    filePathTAS = pathsTas(iName).name;
    filePathTOS = pathsTos(iName).name;  
    [X,Y,VarTAS,Yr,Mth] = loadcmip5out(filePathTAS,'tas');
    [~,~,VarTOS,~,~] = loadcmip5out(filePathTOS,'tos');    
    if sum(ismember(Yr,yrSt:yrEnd)) ~= (yrEnd-yrSt+1)*12 %exclude incomplete runs
        cmipReg(iName).bMat = NaN;
        cmipReg(iName).bIntMat = NaN;
        cmipReg(iName).vPMat = NaN;
        cmipReg(iName).pathTas = filePathTAS;
        cmipReg(iName).pathTos = filePathTOS;
        cmipReg(iName).T = NaN;
        cmipReg(iName).Tp = NaN;
        cmipReg(iName).lags = NaN;
    else
    [~,T] = takespatialavgmonthly(Yr,Mth,Y,VarTAS, 2);
    T = T(Yr >= yrSt & Yr <= yrEnd);
    nino = calcnino34(Yr,Mth,X,Y,VarTOS);
    ninoI = find(f(1).yr == Yr(1));
    predW(ninoI(1):length(nino)+ninoI(1)-1,5) = nino;
    predW = predW(1:length(f(1).yr),:);
    
    
    %Use lag script
    [fLags,~] = searchlags(lMax,predW,yrsI,T,annual);
       %Adjust predictor to reflect chosen lags
    predL = predW;
    for ii = 2:size(predL,2)
        predL(fLags(ii-1)+1:end,ii) = predW(1:end-fLags(ii-1),ii);
    end
    pred = predL(yrsI,:);
    [b,bInt,~,~,stats] = regress(T,pred);
    
    cmipReg(iName).bMat = b;
    cmipReg(iName).bIntMat = bInt;
    cmipReg(iName).vPMat = stats(1);
    cmipReg(iName).pathTas = filePathTAS;
    cmipReg(iName).pathTos = filePathTOS;
    cmipReg(iName).T = T;
    cmipReg(iName).Tp = nansum(repmat(b',[size(pred,1) 1]).* pred,2);
    cmipReg(iName).lags = fLags;
    end
    disp(num2str(cmipReg(iName).lags));
end

%save result
save('/Users/teda/Drive/Research/MATLAB/Solar/paper/mat_files/cmipreghist18822005_20_10_20.mat','cmipReg')
