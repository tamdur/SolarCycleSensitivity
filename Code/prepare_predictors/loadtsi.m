function [mtsi, atsi] = loadtsi
%LOADTSI Store in structure array to specifications used in loadallforcing,
%where dates are from 1850 to 5 years after present. All TSI records
%updated past their end using SOLARIS-HEPPA projections
%
% Ted Amdur
% 9/28/20

mtsi.sources = ["Lean (2000)";...
    "Lean (2000), Wang et al. (2005) background";
    "CMIP6";
    "NRLTSI2"];
atsi.sources = ["Lean (2000)";...
    "Lean (2000), Wang et al. (2005) background";
    "CMIP6";
    "NRLTSI2"];

%Set date vector
cYr = datetime(now,'ConvertFrom','datenum').Year;
aYr = (1850:(cYr+5))';
mYr = repelem(aYr,12);
mth = (1:12)';
mth = repmat(mth,[length(mYr)/12 1]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load all the TSI records
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Incorporate Lean 2000 and WLS 2005
s1gmakewlssolarrecord;
m1 = tsiCycleM;
m1Yr = yrTSIM;
m1Mth = mthTSIM;
m2 = tsiTotM;
m2Yr = yrTSIM;
m2Mth = mthTSIM;
%Incorporate CMIP6 forcing record
m3 = ncread('solarforcing-ref-mon_input4MIPs_solar_CMIP_SOLARIS-HEPPA-3-2_gn_185001-229912.nc','tsi');
m3Yr = ncread('solarforcing-ref-mon_input4MIPs_solar_CMIP_SOLARIS-HEPPA-3-2_gn_185001-229912.nc','calyear');
m3Mth = ncread('solarforcing-ref-mon_input4MIPs_solar_CMIP_SOLARIS-HEPPA-3-2_gn_185001-229912.nc','calmonth');
%Incorporate NRLTSI
m4 = ncread('itsi_ncdc_monthly.nc','TSI');
t4 = ncread('itsi_ncdc_monthly.nc','time');
m4Yr = 1882 + floor(t4./12);
m4Mth = mod(t4,12)+1;

%Prepare array to hold it all
mtsi.tsi = NaN(length(mYr),4);
atsi.tsi = NaN(length(aYr),4);

i1 = find(mYr == m1Yr(1));
mtsi.tsi(i1(1):length(m1)+i1(1)-1,1) = m1;
i2 = find(mYr == m2Yr(1));
mtsi.tsi(i2(1):length(m2)+i2(1)-1,2) = m2;
i3 = find(mYr == m3Yr(1));
mtsi.tsi(i3(1):end,3) = m3(1:length(mYr));
i4 = find(mYr == m4Yr(1));
mtsi.tsi(i4(1):length(m4)+i4(1)-1,4) = m4;

%Update first other records with latest CMIP6 record
mtsi.tsi(i1(1)+length(m1):end,1) = mtsi.tsi(i1(1)+length(m1):end,3)-...
    nanmean(mtsi.tsi(i1(1)+length(m1)-11*12:i1(1)+length(m1),3))+...
    nanmean(mtsi.tsi(i1(1)+length(m1)-11*12:i1(1)+length(m1),1));
mtsi.tsi(i2(1)+length(m2):end,2) = mtsi.tsi(i2(1)+length(m2):end,3)-...
    nanmean(mtsi.tsi(i2(1)+length(m2)-11*12:i2(1)+length(m2),3))+...
    nanmean(mtsi.tsi(i2(1)+length(m2)-11*12:i2(1)+length(m2),2));
mtsi.tsi(i4(1)+length(m4):end,4) = mtsi.tsi(i4(1)+length(m4):end,3)-...
    nanmean(mtsi.tsi(i4(1)+length(m4)-11*12:i4(1)+length(m4),3))+...
    nanmean(mtsi.tsi(i4(1)+length(m4)-11*12:i4(1)+length(m4),4));

%Last, get annual data from above
s = size(mtsi.tsi);
atsi.tsi = squeeze(nanmean(reshape(mtsi.tsi,[12 s(1)/12,s(2)]),1));

%Put dates in
mtsi.yr = mYr;
mtsi.mth = mth;
atsi.yr = aYr;

end



