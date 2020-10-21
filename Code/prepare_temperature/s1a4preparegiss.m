function [X, Y, Yr,Mth, Var] = s1a4preparegiss
%s1a4preparegiss Read in GISS temp for processing
%
% Ted Amdur
% 8/14/20

path = 'gistemp1200_GHCNv4_ERSSTv5.nc';
Y = ncread(path,'lat');
X = ncread(path,'lon');
Var = ncread(path,'tempanomaly');
%Get time variable, convert to year and month
T = ncread(path,'time');
%Note: T is originally in units of days since 1800-1-1, where January
%1850 is entered as 15.5
gissStDate = juliandate(datetime(1800,1,1)); %first date in GISS
jT = T + gissStDate - 1;
dateT = datetime(jT,'convertfrom','juliandate');
Yr = dateT.Year;
Mth = dateT.Month;
end

