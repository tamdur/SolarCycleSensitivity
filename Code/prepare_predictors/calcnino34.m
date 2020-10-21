function [seriesNino, ninoFilt] = calcnino34(Yr,Mth,X,Y,varN)
 %Calculate Nino 3.4 index using recipe from 
% https://climatedataguide.ucar.edu/climate-data/nino-sst-indices-nino-12-3-34-4-oni-and-tni
% Nino X Index computation: (a) Compute area averaged total SST from Niño X region; 
% (b) Compute monthly climatology (e.g., 1950-1979) for area averaged total 
% SST from Niño X region, and subtract climatology from area averaged total 
% SST time series to obtain anomalies; (c) Smooth the anomalies with a 5-month 
% running mean; (d) Normalize the smoothed values by its standard deviation 
% over the climatological period. 
%
% Ted Amdur
% 7/8/2019

%(a) Compute area avg SST for 3.4, which is between 5S-5N, 170-120W
latL = -5; latH = 5; lonL = 360-170; lonH = 360-120;
varN(varN < -1000) = NaN;
varNino = varN(X<=lonH & X >=lonL,Y<=latH & Y>=latL,:);
% seriesNorm = squeeze(nanmean(nanmean(varN(:,Y<=20 & Y >= -20,:),2),1));
seriesNino = squeeze(nanmean(nanmean(varNino,2),1)); %- seriesNorm; %remove warming signal
%(b) Compute monthly climatology for SST, subtract to get anomalies
ninoClim = zeros(12,1);
for iM = 1:12
    mInd = ismember(Mth,iM);
    ninoClim(iM) = nanmean(seriesNino(mInd));
    seriesNino(mInd) = seriesNino(mInd) - ninoClim(iM);
end
%(c) smooth anomalies with a 5-month running mean
seriesNino = movmean(seriesNino,5,'omitnan');
%(d) normalize
seriesNino = seriesNino-nanmean(seriesNino);
seriesNino = seriesNino ./ nanstd(seriesNino);

if nargout > 1
%Perform 8 year high pass to record (from Misios et al. 16, but not using Lanczos which isn't great)
TL = 6*12; flc = 1/TL; %period of low freq cutoff in years (tuned to cut all past 8 years)
Hd1 = designfilt('highpassiir', 'FilterOrder', 15, 'StopBandFrequency', ...
                 flc,'DesignMethod', 'cheby2');
ninoFilt = filtfilt(Hd1,seriesNino);
end
end

