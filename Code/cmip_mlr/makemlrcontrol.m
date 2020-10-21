%Run cmipreg for all cmip control runs, returning set of all estimated 
%sensitivities.
%
% Ted Amdur
% 10/5/20, drawing from generateregcmipmodelshist_1_2020.,

%NOTE: Change cmipreg function pred matrix to use different set of forcings

%Note that recSt and recEnd has to be between 1850,2012 to fit miller
%forcings
clearvars
offset = 3; %Years between samples
yrSt = 1882; yrEnd = 2005;
annual = 0; %Perform annual vs monthly mlr on CMIP runs
tsiChoice = 2;          %1 Lean (2000) 2 WLS (2005) 3 CMIP6 SOLARIS-HEPPA 4. NRLTSI2
detrendTSI = 1;         %1 to detrend TSI over interval, 0 otherwise
anthropogenicChoice = 1;%1 Miller et al. (2014) 2 Dessler Forster (2018) 3 RCP 6.0 effCO2
volcanismChoice = 1;    %1 Sato 550nm mean aod 2. CMIP6 550 nm mean aod

%hist files are stored in: = dir('/Users/teda/Data/Odyssey/tamdur/*.nc');
pathsTas = dir('/Users/teda/Data/Odyssey/tascontrol/*.nc');
pathsTos = dir('/Users/teda/Data/Odyssey/toscontrol/*.nc');

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
    %[~,~,VarTAS,~,~] = loadcmip5out(filePathTAS,'tas');
    if iName == 17
        VarTAS = VarTAS(:,:,1:length(Yr));
    end
    if sum(ismember(Yr,Yr(1):Yr(end))) ~= (Yr(end)-Yr(1)+1)*12 %exclude incomplete runs
        cmipReg(iName).bMat = NaN;
        cmipReg(iName).vPMat = NaN;
        cmipReg(iName).pathTas = filePathTAS;
        cmipReg(iName).pathTos = NaN;
        cmipReg(iName).T = NaN;
        cmipReg(iName).Tp = NaN;
        cmipReg(iName).lags = NaN;
        cmipReg(iName).numSub = 0;
    else
    [~,TAll] = takespatialavgmonthly(Yr,Mth,Y,VarTAS, 2);
    ninoAll = calcnino34(Yr,Mth,X,Y,VarTAS);
    ninoI = find(f(1).yr == yrSt);
    ninoRep = ninoI(1):length(yrsI)+ninoI(1)-1;
    
    %Make blocks within this long time series
    ind = 1;
    numSub = length(Yr(1)+1:offset:Yr(end) - (yrEnd-yrSt));
    bMat = NaN(numSub,5); vPMat = NaN(numSub,1); lagMat = NaN(numSub,4);
    for iBlock = Yr(1)+1:offset:Yr(end) - (yrEnd-yrSt)
        predL = predW;
        iYrs = find(Yr >= iBlock & Yr <= iBlock + (yrEnd-yrSt));
        T = TAll(iYrs);
        nino = ninoAll(iYrs);
        predL(ninoRep,5) = nino;
        %Use lag script
        [fLags,~] = searchlags(lMax,predL,yrsI,T,annual);
           %Adjust predictor to reflect chosen lags
        for ii = 2:size(predL,2)
            predL(fLags(ii-1)+1:end,ii) = predL(1:end-fLags(ii-1),ii);
        end
        pred = predL(yrsI,:);
        [b,bInt,~,~,stats] = regress(T,pred);
        bMat(ind,:) = b;
        vPMat(ind) = stats(1);
        lagMat(ind,:) = fLags;
        ind = ind + 1;
    end
    
    cmipReg(iName).bMat = bMat;
    cmipReg(iName).vPMat = vPMat;
    cmipReg(iName).pathTas = filePathTAS;
    cmipReg(iName).pathTos = NaN;
    cmipReg(iName).T = TAll;
    cmipReg(iName).lags = lagMat;
    cmipReg(iName).numSub = numSub;
    end
    disp(['Run ' num2str(iName) ' of ' num2str(length(pathsTas)) ' completed'])
end

%save result
save('/Users/teda/Drive/Research/MATLAB/Solar/paper/mat_files/cmipmlrcontrol18822005_16_10_20.mat','cmipReg')
