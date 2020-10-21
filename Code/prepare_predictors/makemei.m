function [mei,yrMei,mthMei] = makemei
%MAKEMEI Make MEI record for loadenso.m

%Generate MEI_ext variable
filename = '../Data/ENSO/mei_ext.txt';
delimiter = '\t';
startRow = 8;
formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';
fileID = fopen(filename,'r');
textscan(fileID, '%[^\n\r]', startRow-1, 'WhiteSpace', '', 'ReturnOnError', false,...
    'EndOfLine', '\r\n');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', ...
    'string', 'ReturnOnError', false);
fclose(fileID);
%Make array with extended MEI
%Format: years from 1871-2005, values as standard deviations from mean
%columns: YEAR    DECJAN  JANFEB  FEBMAR  MARAPR  APRMAY  MAYJUN  JUNJUL  JULAUG ...
%AUGSEP  SEPOCT  OCTNOV  NOVDEC
meiExt = [dataArray{1:end-1}];
meiOld = meiExt(:,2:end)';
meiOld = meiOld(:);
meiYr = (1871:2005)'; 
meiYr = repelem(meiYr,12);
meiMth = (1:12)';
meiMth = repmat(meiMth,[length(1871:2005) 1]);
clearvars filename delimiter startRow formatSpec fileID dataArray ans;

%Load new MEI
ext = ncread('imeiv2.nc','MEI');
extT = ncread('imeiv2.nc','time');
mthExt = mod(extT,12)+1;
yrExt = 1979 + floor(extT./12);




%Stitch together the extended MEI with MEI to get the ENSO predictor
%Uses MEI 1950-2018, mei Ext earlier
cYr = datetime(now,'ConvertFrom','datenum').Year;
aYr = (1850:(cYr+5))';
mYr = repelem(aYr,12);
mth = (1:12)';
mth = repmat(mth,[length(mYr)/12 1]);
mei = NaN(length(mYr),1);
iOld = find(mYr==1871);
mei(iOld(1):iOld(1)+length(meiOld)-1) = meiOld;
iExt = find(mYr == 1979);
mei(iExt(1):iExt(1)+length(ext)-1) = ext;
yrMei = mYr; mthMei = mth;
end

