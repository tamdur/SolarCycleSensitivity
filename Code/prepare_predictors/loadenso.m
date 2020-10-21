function [me, ae] = loadenso
%LOADENSO Store in structure array with commonly indexed years/months
%different valid predictors for ENSO. Note that monthly records are
%averaged to create an alternative annual record, but the converse is not
%performed. 
%
% Ted Amdur
% 9/29/20
me.sources = ["Kaplan 3.4";...
    "Kaplan 3.4 8yrHP";...
    "MEI"; ...
    "ESRL 3.4"];
ae.sources = [me.sources; "MEI prior DJF"];


%Set date vector
cYr = datetime(now,'ConvertFrom','datenum').Year;
aYr = (1850:(cYr+5))';
mYr = repelem(aYr,12);
mth = (1:12)';
mth = repmat(mth,[length(mYr)/12 1]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load all the ENSO records
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Make ENSO using Kaplan SST
xK = ncread('sst.mon.anom.nc','lon');
yK = ncread('sst.mon.anom.nc','lat');
tK = ncread('sst.mon.anom.nc','time');
varK = ncread('sst.mon.anom.nc','sst');
tKStart = juliandate(datetime(1800,1,1));
jT = tK + tKStart;
dateT = datetime(jT,'convertfrom','juliandate');
m1Yr = dateT.Year; m2Yr = m1Yr;
MthK = dateT.Month;
[m1, m2] = calcnino34(m1Yr,MthK,double(xK),double(yK),double(varK));

%Load MEI
m3 = makemei; %Already in 1850 to 5 years past present formatting

%Load ERSST
m4 = ncread('iersst_nino3.4a_rel.nc','Nino3.4r');
m4T = ncread('iersst_nino3.4a_rel.nc','time');
m4Mth = mod(m4T,12)+1;
m4Yr = 1854 + floor(m4T./12);

%Prepare array to hold it all
me.enso = NaN(length(mYr),4);
ae.enso = NaN(length(aYr),5);

i1 = find(mYr == m1Yr(1));
me.enso(i1(1):length(m1)+i1(1)-1,1) = m1;
i2 = find(mYr == m2Yr(1));
me.enso(i2(1):length(m2)+i2(1)-1,2) = m2;
me.enso(:,3) = m3;
i4 = find(mYr == m4Yr(1));
me.enso(i4(1):length(m4)+i4(1)-1,4) = m4;

%Last, get annual data from above
s = size(me.enso);
ae.enso(:,1:4) = squeeze(nanmean(reshape(me.enso,[12 s(1)/12,s(2)]),1));

%Put dates in
me.yr = mYr;
me.mth = mth;
ae.yr = aYr;

%Last, average for DJF
ind = 1;
djfInd = find(ismember(me.mth,[12 1 2]));
for ii = me.yr(13):me.yr(end)
    thisYr = find(me.yr==ii);
    thisYr = [thisYr(1)-1; thisYr(1:end-1)];
    ensoI = intersect(djfInd,thisYr);
    aenso(ind) = nanmean(me.enso(ensoI,3),1);
    ind = ind + 1;
end
ae.enso(2:end,5) = aenso;



end



