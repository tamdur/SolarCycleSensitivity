function [ma, aa] = loadanth
%LOADALLANTH Store in structure array with commonly indexed years/months
%different valid predictors for anthropological forcing.
%performed. 
%
% Ted Amdur
% 9/5/19

ma.sources = ["Miller et al. (2014)"; 
    "Dessler Forster (2018)";
    "RCP effective co2 concentration"];
aa.sources = ma.sources;

%Set date vector
cYr = datetime(now,'ConvertFrom','datenum').Year;
aYr = (1850:(cYr+5))';
mYr = repelem(aYr,12);
mth = (1:12)';
mth = repmat(mth,[length(mYr)/12 1]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load all the Anthropogenic records
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Load Miller Anthro forcing
filename = '../Data/Anthropogenic/Fe_H11_1880-2011.txt';
startRow = 14;
formatSpec = '%4f%8f%8f%8f%8f%8f%8f%8f%8f%8f%f%[^\n\r]';
fileID = fopen(filename,'r');
textscan(fileID, '%[^\n\r]', startRow-1, 'WhiteSpace', '', 'ReturnOnError', ...
    false, 'EndOfLine', '\r\n');
dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', ...
    'TextType', 'string', 'ReturnOnError', false);
fclose(fileID);
a1Yr = dataArray{:, 1};
WMGHGs = dataArray{:, 2};
O3 = dataArray{:, 3};
StrH2O = dataArray{:, 4};
ReflAer = dataArray{:, 5};
AIE = dataArray{:, 6};
BC = dataArray{:, 7};
SnowAlb = dataArray{:, 8};
StrAer = dataArray{:, 9};
Solar = dataArray{:, 10};
LandUse = dataArray{:, 11};
a1 = WMGHGs+O3+StrH2O+ReflAer+AIE+BC+SnowAlb+LandUse;

a2 = ncread('iAnthropogenic_total_ERF.nc','Anthropogenic_total_ERF');
a2Yr = ncread('iAnthropogenic_total_ERF.nc','time');
a2Yr = a2Yr + 1750;

%Load ghg data from CMIP5 effective co2 concentration
% mthGHG = ncread('cmip5_GhgCon_rcp45_1850-2020_interpol_mm.nc','time');
% co2 = ncread('cmip5_GhgCon_rcp45_1850-2020_interpol_mm.nc','co2');
% %Note: the above is the effective CO2 interpolated concentration from all
% %GHGs
% ghgStDate = datetime(1765,6,1); %first date in record
filename = '../Data/Anthropogenic/iRCP6_CO2EQ.dat';
startRow = 6;
formatSpec = '%4f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '',...
    'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false,...
    'EndOfLine', '\r\n');
fclose(fileID);
a3Yr = dataArray{:, 1};
a3 = dataArray{:, 2};
%convert to radiative forcing, following prescription given in paper
a3 = 3.71/log(2)*log(a3/278);

%Prepare array to hold it all
aa.anth = NaN(length(aYr),3);
ma.anth = NaN(length(mYr),3);

i1 = find(aYr == a1Yr(1));
aa.anth(i1(1):length(a1)+i1(1)-1,1) = a1;
aa.anth(1:i1(1)-1,1) = 0;
i2 = find(a2Yr == aYr(1));
aa.anth(1:length(a2(i2:end)),2) = a2(i2:end);
i3 = find(a3Yr == aYr(1));
aa.anth(:,3) = a3(i3:length(aYr)+i3-1);

%Update records with the RCP 6.0 scencario up through end, using past ten
%years to normalize
aa.anth(i1(1)+length(a1):end,1) = aa.anth(i1(1)+length(a1):end,3)-...
    nanmean(aa.anth(i1(1)+length(a1)-10:i1(1)+length(a1),3))+...
    nanmean(aa.anth(i1(1)+length(a1)-10:i1(1)+length(a1),1));
aa.anth(length(a2(i2:end)):end,2) = aa.anth(length(a2(i2:end)):end,3)-...
    nanmean(aa.anth(length(a2(i2:end))-10:length(a2(i2:end)),3))+...
    nanmean(aa.anth(length(a2(i2:end))-10:length(a2(i2:end)),2));

%Last, get monthly data from above
ma.anth = repelem(aa.anth,12,1);

%Put dates in
ma.yr = mYr;
ma.mth = mth;
aa.yr = aYr;



end



