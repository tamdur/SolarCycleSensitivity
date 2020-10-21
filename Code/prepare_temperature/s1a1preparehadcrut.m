function [X, Y, Yr,Mth, Var] = s1a1preparehadcrut(varargin)
%First step of reading in hadcrut data for processing. From camptung_hadcrut.m
%Ted Amdur
%3/5/2019
%
%   Outputs:
%           X = longitude of gridpoint
%           Y = latitude of gridpoint
%           Yr = Calendar Year corresponding to ith slice of Var
%           Mth = Month (1-12) corresponding to ith slice of Var
%           Var = Reanalysis record, of dimensions [X,Y,T] 

%Load data from HadCRUT
if nargin < 1
    path = 'HadCRUT.4.6.0.0.median.nc';
    tVar = 'temperature_anomaly';
    timeVar = 'time';
else %For use with HadCRUT3
    path = char(varargin);
    tVar = 'temp';
    timeVar = 't';
end
Y = ncread(path,'latitude');
X = ncread(path,'longitude');
Var = ncread(path,tVar);
%Get time variable, convert to year and month
T = ncread(path,timeVar);
%Note: T is originally in units of days since 1850-1-1, where January
%1850 is entered as 15.5
hadStDate = juliandate(datetime(1850,1,1)); %first date in HadCRUT
jT = T + hadStDate - 1;
dateT = datetime(jT,'convertfrom','juliandate');
Yr = dateT.Year;
Mth = dateT.Month;

end