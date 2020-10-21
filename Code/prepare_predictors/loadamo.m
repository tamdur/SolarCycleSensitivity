function [mamo,aamo] = loadamo
%LOADAMO Return structure array with commonly indexed years/months
%the PSL AMO index. Note that annual record generated by averaging monthly 
%monthly records. Indexed beginning in 1850 in accordance with 
%loadallforcing formatting.
%
% Ted Amdur
% 9/29/20
% Metadata from file:
%   AMO unsmoothed, detrended from the Kaplan SST V2
%   Calculated at NOAA PSL1
%   http://www.psl.noaa.gov/data/timeseries/AMO/

%Read in AMO
%note that as of this comment, ENSO 3.4 is stored as 'nino34.anom.data'
aamo.sources = ["HadISST";"ERSSTv5"];
mamo.sources = aamo.sources;


%Set date vector
cYr = datetime(now,'ConvertFrom','datenum').Year;
aYr = (1850:(cYr+5))';
mYr = repelem(aYr,12);
mth = (1:12)';
mth = repmat(mth,[length(mYr)/12 1]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load all the AMO records
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m1 = ncread('iamo_hadsst_ts.nc','AMO');
mT = ncread('iamo_hadsst_ts.nc','time');
m1Yr = 1853 + floor(mT./12);
m2 = ncread('iamo_ersst.nc','AMO');
mT = ncread('iamo_ersst.nc','time');
m2Yr = 1880 + floor(mT./12);

%Prepare array to hold it all
mamo.amo = NaN(length(mYr),2);
aamo.amo = NaN(length(aYr),2);

i1 = find(mYr == m1Yr(1));
mamo.amo(i1(1):length(m1)+i1(1)-1,1) = m1;
i2 = find(mYr == m2Yr(1));
mamo.amo(i2(1):length(m2)+i2(1)-1,2) = m2;

%Last, get annual data from above
s = size(mamo.amo);
aamo.amo = (squeeze(nanmean(reshape(mamo.amo,[12 s(1)/12,s(2)]),1)));

%Put dates in
mamo.yr = mYr;
mamo.mth = mth;
aamo.yr = aYr;

end

