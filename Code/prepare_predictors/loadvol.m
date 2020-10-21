function [mv, av] = loadvol
%LOADALLVOL Store in structure array with commonly indexed years/months
%different valid predictors for volcanism. Note that monthly records are
%averaged to create an alternative annual record, but the converse is not
%performed. 
%
% Ted Amdur
% 9/5/19


mv.sources = ["Sato  550 nm mean aod";"CMIP6 550 nm mean aod"];
av.sources = ["Sato  550 nm mean aod";"CMIP6 550 nm mean aod"];


%Set date vector
cYr = datetime(now,'ConvertFrom','datenum').Year;
aYr = (1850:(cYr+5))';
mYr = repelem(aYr,12);
mth = (1:12)';
mth = repmat(mth,[length(mYr)/12 1]);


% From https://data.giss.nasa.gov/modelforce/strataer/
filename = '/Users/teda/Drive/Research/Data/Forcing/Lean_Rind/tau.line_2012.12.txt';
startRow = 5;
formatSpec = '%8f%8f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '',...
    'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false,...
    'EndOfLine', '\r\n');
fclose(fileID);
v1YrF = dataArray{:, 1};
v1Yr = floor(v1YrF);
v1 = dataArray{:,2};



%read in relevant data
alt = ncread('CMIP_1850_2014_extinction_550nm_strat_only_v3.nc','altitude');
lat = ncread('CMIP_1850_2014_extinction_550nm_strat_only_v3.nc','latitude');
month = ncread('CMIP_1850_2014_extinction_550nm_strat_only_v3.nc','month');
ext550 = ncread('CMIP_1850_2014_extinction_550nm_strat_only_v3.nc','ext550');

%First, convert into year and month variables
%Note: T is originally in units of months since 1850-1-1, where January
%1850 is entered as 1. Note I believe the metadata from the nc file is
%incorrect, saying that it is months from jan 1960.
v2Yr = floor((month-1)/12) + 1850; %Return year
v2Mth = mod(month-1,12)+1; %Add 1 such at January is month 1 (vs month 0)


%Take the weighted latitudinal average
wVec = cos((pi/180).*lat)';
wVec = wVec ./ (sum(wVec)); %Make it add up to 1
for iL = 1:36
    ext550(:,:,iL) = ext550(:,:,iL).*wVec(iL);
end
ext550 = nansum(ext550,3); %Take lat average

%Next, get the AOD timeseries by integrating the extinction coefficient
%above 15km
stratLvl = alt>=15;
v2 = nansum(ext550(:,stratLvl),2);

%Prepare array to hold it all
av.vol = NaN(length(aYr),2);
mv.vol = NaN(length(mYr),2);

mv.vol(1:length(v1),1) = v1;
mv.vol(1:length(v2),2) = v2;


%Update through end using last 10 years of data
mv.vol(length(v1)+1:end,1) = nanmean(mv.vol(length(v1)-10*12:length(v1),1));
mv.vol(length(v2)+1:end,2) = nanmean(mv.vol(length(v2)-10*12:length(v2),2));

%Convert to units of radiative forcing
delF = -26; %Conversion factor for 550nm AOD as determined by Schmidt et al. 18,
            %in line with Hansen et al 05 estimate of -25;
mv.vol = mv.vol .* delF;

%Last, get annual data from above
s = size(mv.vol);
av.vol = squeeze(nanmean(reshape(mv.vol,[12 s(1)/12,s(2)]),1));

%Put dates in
mv.yr = mYr;
mv.mth = mth;
av.yr = aYr;

end

