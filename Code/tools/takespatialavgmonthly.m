function [year, series] = takespatialavgmonthly(yr,mth,Y,varN, method)
%Version of takespatialavg that produces a monthly timeseries. For use with NCEP,
%CMIP output, Berkeley spatial temperature dataset, HadCRUT. Factors in missing
%values with spatial weighting, subtracts away monthly climatology
%
% Ted Amdur
% 6/15/19
%
%   INPUTS: 
%           method: (1) dumb, simple average. (2) latitude weighted

%Weighted average for each band is given by the sum of sum of lat band anomalies 
%times the spatial weight. Total Result is then divided by area of non-nan
%cells

%Generate weighting vector for latitudes (assuming Y in degrees)
if method == 1
    wVec = ones(length(Y),1)';
else
    wVec = cos((pi/180).*Y)';
end
wVec = wVec ./ (sum(wVec)); %Make it add up to 1

series = zeros(size(mth,1),1);
for ii = 1:length(yr)
    %Select a given month
    varY = varN(:,:,ii);
    %Sum anomalies in each band
    latSum = nansum(varY,1);
    %weigh by the latitude
    latSum = latSum.*wVec;
    %get total area of non-nan cells
    obs = double(~isnan(varY));
    weighObs = sum(sum(obs.*wVec));
    %divide by area of valid cells
    yT = nansum(latSum)/weighObs;
    series(ii) = yT;
    if ii == 1400
        a = 1;
    end
end
year = yr; 
% series = series';

seriesClim = zeros(12,1);
for iM = 1:12
    mInd = ismember(mth,iM);
    seriesClim(iM) = nanmean(series(mInd));
    series(mInd) = series(mInd) - seriesClim(iM);
end

%Make sure series is oriented in correct (NX1 vector) direction
if size(series,2) > 1
    series = series';
end