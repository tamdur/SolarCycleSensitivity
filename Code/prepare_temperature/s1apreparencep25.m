function [X, Y, Yr,Mth, Var] = s1apreparencep25
%First step of reading in NCEP data for processing. from schprocess_NCEP.m
%Ted Amdur
%2/12/2019
%   Inputs: 
%       fileName: String containing the file location of NCEP file. Current
%       file locations are 'NCEP_NH_sfcpressure.nc', 'NCEP_Tsfc.nc'
%       'NCEP_NH_accumsnow.nc', 'NCEP_NH_sfctemp.nc', and 'NCEP_NH_cld.nc'
%       type: String containing the type of data being unpacked. Options:
%           -snow: snow depth (surface)
%           -pressure: surface pressure
%           -cloud: total cloud fraction
%           -temp: surface temperature
%
%   Outputs:
%           X = longitude of gridpoint
%           Y = latitude of gridpoint
%           Yr = Calendar Year corresponding to ith slice of Var
%           Mth = Month (1-12) corresponding to ith slice of Var
%           Var = Reanalysis record, of dimensions [X,Y,T] 

fileName = 'Tsfc25.nc';
Var = ncread(fileName,'air');
Var = squeeze(Var); %squeeze 1-D z field 

X = ncread(fileName,'lon');
Y = ncread(fileName,'lat');

%Get time variable, convert to year and month
T = ncread(fileName, 'time');
%Note: T is originally in units of hours since since 1800-1-1
T = T /24; %make it into days since 1800
stDay = juliandate(datetime(1800,1,1,0,0,0));
T = T + stDay;
time = datetime(T,'ConvertFrom','juliandate');
Yr = time.Year; %Return year
Mth = time.Month; %Add months to make all entries positive

end