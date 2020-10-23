function cmipU = cmip_T_uncertainty
%calculate the temperature uncertainty due to internal variability by
%finding the mismatch between CMIP5 historical outputs and predicted
%temperature
% Ted Amdur
% 2/10/20


%First, load the file of predicted temperatures across historical runs
load likeshist_2_10_20.mat;
histLik = histLik([1:35 37:end]); %Address broken run

%Select the runs for which the predictors match the actual CMIP5 forcings.
%These are:
%solar: WLS 2, Anthro: Miller 1, vol: Sato 1, ENSO: internal 3.4 1
%Then, for each extract the temperature time series and predicted
%temperature time seriesc
cmipU.T = [];
cmipU.TP = [];
cmipU.dif = [];
for ii = 1:length(histLik)
    hI = histLik(ii).predM;
    selI = find(hI(:,1) == 2 & hI(:,2) == 1 & hI(:,3) == 1 & ...
        hI(:,4) == 1);
    T = histLik(ii).T;
    TP = histLik(ii).tpMat(selI,:);
    cmipU.T = [cmipU.T; T'];
    cmipU.TP = [cmipU.TP; TP];
    cmipU.dif = [cmipU.dif; (T-TP')'];
end
cmipU.SD = std(cmipU.dif,0,1);
end
    