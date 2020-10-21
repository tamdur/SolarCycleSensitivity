%Process the AOD files from CMIP6 (located with readme in
%Drive/Data/Volcanic) into global mean 550nm aod and associated month, year
%
% Ted Amdur
% 10/22/19

%read in relevant data
alt = ncread('CMIP_1850_2014_extinction_550nm_strat_only_v3.nc','altitude');
lat = ncread('CMIP_1850_2014_extinction_550nm_strat_only_v3.nc','latitude');
month = ncread('CMIP_1850_2014_extinction_550nm_strat_only_v3.nc','month');
ext550 = ncread('CMIP_1850_2014_extinction_550nm_strat_only_v3.nc','ext550');

%First, convert into year and month variables
%Note: T is originally in units of months since 1850-1-1, where January
%1850 is entered as 1. Note I believe the metadata from the nc file is
%incorrect, saying that it is months from jan 1960.
Yr = floor((month-1)/12) + 1850; %Return year
Mth = mod(month-1,12)+1; %Add 1 such at January is month 1 (vs month 0)


%Take the weighted latitudinal average
wVec = cos((pi/180).*lat)';
wVec = wVec ./ (sum(wVec)); %Make it add up to 1
for iL = 1:36
    ext550(:,:,iL) = ext550(:,:,iL).*wVec(iL);
end
ext550 = nansum(ext550,3); %Take lat average

%Next, get the AOD timeseries by integrating the extinction coefficient
%above 15km
stratLvl = find(alt>=15);
extS = nansum(ext550(:,stratLvl),2);

% %Save as a simple structure over 1850 to end of 2018
% f = s1hmakeforcings;
% cmipvol.vol = NaN(length(f.mYr),1);
% cmipI = find(ismember(f.mYr,Yr));
% cmipvol.vol(cmipI) = extS(ismember(Yr,f.mYr));
% cmipvol.mth = f.mMth; cmipvol.yr = f.mYr;
% save('cmip6vol.mat','cmipvol')


