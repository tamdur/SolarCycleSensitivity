%Generate and save the WLS 2005 solar forcing record used in CMIP5 model
%studies and incorporated by all studies shown in paper, with exception of
%reconstructions for paleo analysis (as of 7/31/19)



%-------------------------------------------------------------------
% Monthly TSI
%-------------------------------------------------------------------
filename = '../Data/TSI/TSI_WLS_mon_1882_2008.txt';
startRow = 4;
formatSpec = '%8f%8f%15f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '',...
    'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false,...
    'EndOfLine', '\r\n');

fclose(fileID);

yrTSIM = dataArray{:, 1};
mthTSIM = dataArray{:, 2};
yrFTSIM = yrTSIM + (mthTSIM-1)./12;
tsiCycleM = dataArray{:, 3};
tsiTotM = dataArray{:, 4};
clearvars filename startRow formatSpec fileID dataArray ans;

%-------------------------------------------------------------------
% Yearly TSI
%-------------------------------------------------------------------
filename = '../Data/TSI/TSI_WLS_ann_1610_2008.txt';
startRow = 4;
formatSpec = '%8f%15f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', ...
    'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false,...
    'EndOfLine', '\r\n');
fclose(fileID);
yrTSIA = dataArray{:, 1};
tsiCycleA = dataArray{:, 2};
tsiTotA = dataArray{:, 3};

clearvars filename startRow formatSpec fileID dataArray ans;



