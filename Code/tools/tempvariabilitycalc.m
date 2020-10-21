function stdSeries = tempvariabilitycalc(stYr,endYr,blkS)
% Use HadCRUT4 ensemble to calculate the uncertainty in temperature
% Ted Amdur
% 1/25/20
[yrE,ens] = loadensemble;
ens = ens(yrE >= stYr & yrE <= endYr,:);
%Get the empirical standard deviation 
subN = floor(size(ens,1)/blkS);
stdSeries = NaN(subN,1);
for ii = 1:subN
    %Calculate the standard deviation of the block
    interval = 1+(ii-1)*blkS:ii*blkS;
    ensBlk = mean(ens(interval,:),1);
    stdSeries(ii) = std(ensBlk,0,2);
end
end


function [yrE,ens] = loadensemble
%Return a vector of years, NX100 matrix holding 100 ensemble timeseries
paths = dir('../Data/Temperature/Time_series/HadCRUT.4.6.0.0.monthly_ns_avg_realisations/');
paths = paths(3:end);
filePath = paths(1).name;
filePath = ['../Data/Temperature/Time_series/HadCRUT.4.6.0.0.monthly_ns_avg_realisations/' filePath];
hadCRUT = s1a1preparehadcrutmonthly(filePath);
date = hadCRUT{:,1};
yrE = date.Year;
ens = NaN(length(yrE),100);
ens(:,1) = hadCRUT{:,2};
for iName = 2:length(paths)
    filePath = paths(iName).name;
    filePath = ['../Data/Temperature/Time_series/HadCRUT.4.6.0.0.monthly_ns_avg_realisations/' filePath];
    hadCRUT = s1a1preparehadcrutmonthly(filePath);
    ens(:,iName) = hadCRUT{:,2};
end
end