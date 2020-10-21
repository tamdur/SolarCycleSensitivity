function [yrSeries,series] = s1a3prepareberkeley%Read in Berkeley time series for processing.%Ted Amdur%9/30/20 using import script from matlabfilename = '/Users/teda/Drive/Research/MATLAB/JCLI_Sun_tools/Data/Temperature/Time_series/Berkeley_Land_and_Ocean_complete.txt';startRow = 78;cYr = datetime(now,'ConvertFrom','datenum').Year;endRow = (cYr-1850).* 12 + 78-1; %Set to reflect the year this is downloadedformatSpec = '%1C%5f%6f%10f%[^\n\r]';fileID = fopen(filename,'r');dataArray = textscan(fileID, formatSpec, endRow-startRow+1, 'Delimiter', '', ...    'WhiteSpace', '', 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines', startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');fclose(fileID);ber = table(dataArray{1:end-1}, 'VariableNames', ...    {'VarName1','Year','Month','TAnomaly'});%% Clear temporary variablesclearvars filename startRow endRow formatSpec fileID dataArray ans;yrSeries = ber.Year;series = ber.TAnomaly;end