function [mipo,aipo] = loadipo
%LOADIPO Store in structure array with commonly indexed years/months
%the PSL PDO index. Note that monthly records are
%averaged to create an annual record.
%
% Ted Amdur
% 8/14/20
%
% Metadata from file:
% TPI unfiltered (IPO Tripole Index) of Henley et al. (2015)
% HadISST 1.1
% created at NOAA PSL
% If you use this data please cite: 
% Henley, B.J., Gergis, J., Karoly, D.J., Power, S.B., Kennedy, J., 
% & Folland, C.K. (2015). A Tripole Index for the Interdecadal Pacific Oscillation. 
%Climate Dynamics, 45(11-12), 3077-3090. doi:10.1007/s00382-015-2525-1
% See https://psl.noaa.gov/data/timeseries/IPOTPI/
% Updated date at NOAA/PSL
% Produced Aug 13 2020
% units="degC"

%Read in IPO
fileName = 'ipo.tpi.timeseries.hadisst11.data';

aipo.sources = "HadISST";
mipo.sources = aipo.sources;

formatSpec = '%5f%9f%8f%8f%8f%8f%8f%8f%8f%8f%8f%8f%f%[^\n\r]';

%% Open the text file.
fileID = fopen(fileName,'r');

%% Read columns of data according to the format.
firstLine = textscan(fileID, formatSpec, 1,'Delimiter', '', 'WhiteSpace', ...
    '', 'TextType', 'string',  'ReturnOnError', false);
lenRead = 2018-firstLine{1}+1; %Other predictors go through 2018
dataArray = textscan(fileID, formatSpec, lenRead,'Delimiter', '', ...
    'WhiteSpace', '', 'TextType', 'string',  'ReturnOnError', false);

%% Close the text file.
fclose(fileID);

%% Create output variable
ipoM = [dataArray{1:end-1}];
ipoM(ipoM < -10) = NaN;
%% Clear temporary variables
clearvars formatSpec fileID dataArray ans;

m1 = ipoM(:,2:end)';
m1 = m1(:);
m1Yr = ipoM(:,1); 
m1Yr = repelem(m1Yr,12);
m1Mth = (1:12)';
m1Mth = repmat(m1Mth,[size(ipoM,1) 1]);

%Set date vector
cYr = datetime(now,'ConvertFrom','datenum').Year;
aYr = (1850:(cYr+5))';
mYr = repelem(aYr,12);
mth = (1:12)';
mth = repmat(mth,[length(mYr)/12 1]);
%Prepare array to hold it all
mipo.ipo = NaN(length(mYr),1);
aipo.ipo = NaN(length(aYr),1);

i1 = find(mYr == m1Yr(1));
mipo.ipo(i1(1):length(m1)+i1(1)-1,1) = m1;

%Last, get annual data from above
s = size(mipo.ipo);
aipo.ipo = (squeeze(nanmean(reshape(mipo.ipo,[12 s(1)/12,s(2)]),1)))';

%Put dates in
mipo.yr = mYr;
mipo.mth = mth;
aipo.yr = aYr;


end

