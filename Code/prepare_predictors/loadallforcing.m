function f= loadallforcing(annual)
%Create forcings to be used for any MLR analysis 
%
%   OUTPUTS:
%            f = structure array. Each iteration of f reflects a category of 
%               forcing with following fields:
%               yr: year vector corresponding to fields. Range up to 5
%               years from 1850 to past call of structure
%               mth: Only if annual = 0. Month vector corresponding to the month of each
%               entry. Range of 1-12 spanning 1850 to 5 years past call of
%               structure
%               yrF: Only if annual = 0. fractional year over range of yr
%               preds: pred array of size NXM, where N is number of
%               observations corresponding to length of yr vector, and M is
%               number of predictors. Units of radiative forcing with
%               exception of internal indices such as ENSO and AMO
%               type: type of predictor contained within, char
%               sources: sources of M predictors, char
%               NOTE: ADD OTHER FORCINGS
%
% Ted Amdur
% 9/28/20

%First, make time vecs
cYr = datetime(now,'ConvertFrom','datenum').Year;

yr = (1850:(cYr+5))';
if ~annual
    yr = repelem(yr,12);
    mth = (1:12)';
    mth = repmat(mth,[length(yr)/12 1]);
    yrF = yr+ (mth-1)/12;
end
    
[tsim,tsia] = loadtsi;
[ghgm,ghga] = loadanth;
[volm,vola] = loadvol;
[ensom,ensoa] = loadenso;
[amom,amoa] = loadamo;
[ipom,ipoa] = loadipo;

for ii = 1:7
    f(ii).yr = yr;
    if ~annual
        f(ii).mth = mth;
        f(ii).yrF = yrF;
    end
end
f(1).type = 'tsi';
if annual
    f(1).preds = tsia.tsi;
    f(1).sources = tsia.sources;
else
    f(1).preds = tsim.tsi;
    f(1).sources = tsim.sources;
end

f(2).type = 'anthropogenic';
if annual
    f(2).preds = ghga.anth;
    f(2).sources = ghga.sources;
else
    f(2).preds = ghgm.anth;
    f(2).sources = ghgm.sources;
end

f(3).type = 'vol';
if annual
    f(3).preds = vola.vol;
    f(3).sources = vola.sources;
else
    f(3).preds = volm.vol;
    f(3).sources = volm.sources;
end

f(4).type = 'enso';
if annual
    f(4).preds = ensoa.enso;
    f(4).sources = ensoa.sources;
else
    f(4).preds = ensom.enso;
    f(4).sources = ensom.sources;
end

f(5).type = 'amo';
if annual
    f(5).preds = amoa.amo;
    f(5).sources = amoa.sources;
else
    f(5).preds = amom.amo;
    f(5).sources = amom.sources;
end

f(6).type = 'ipo';
if annual
    f(6).preds = ipoa.ipo;
    f(6).sources = ipoa.sources;
else
    f(6).preds = ipom.ipo;
    f(6).sources = ipom.sources;
end

f(7).type = 'other';
% if annual
%     f(7).preds = otha.other;
%     f(7).sources = otha.sources;
% else
%     f(7).preds = othm.other;
%     f(7).sources = othm.sources;
% end

end


