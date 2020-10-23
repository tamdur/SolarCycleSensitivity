%Draw from EnergyBand in ph_scripts folder to estimate the amount of
%variability in the decadal variability band from each predictor.
%
% Ted Amdur
% 4/14/20, update 10/22/20
stYr = 1882;
endYr = 2019;
load loadallforcingstruct.mat
f = fMon;
% [tsi,~] = loadalltsi;
% [ghg,~] = loadallanth(1);
% [vol,~] = loadallvol;
% [enso,~] = loadallenso;

%Make pred for ASH20
tsiSub = f(1).preds(f(1).yr >= stYr & f(1).yr <= endYr,3);
tsiSub = detrend(tsiSub,'linear');
ghgSub = f(2).preds(f(2).yr >= stYr & f(2).yr <= endYr,2);
volSub = f(3).preds(f(3).yr >= stYr & f(3).yr <= endYr,2);
ensoSub = f(4).preds(f(4).yr >= stYr & f(4).yr <= endYr,3);


dt=1/12;


[freqS,PS]=fftPH(tsiSub,dt,[],0);
[freqA,PA]=fftPH(ghgSub,dt,[],0);
[freqV,PV]=fftPH(volSub,dt,[],0);
[freqE,PE]=fftPH(ensoSub,dt,[],0);

plS=find(freqS>=1/14 & freqS<=1/8);  

figure(1); clf; hold on;
plot(freqS,PS); logPH;
plot(freqS(plS),PS(plS),'r');
hold on
plot(freqA,PA,'LineStyle',':');
plot(freqA(plS),PA(plS),'r','LineStyle',':');
hold on
plot(freqV,PV,'LineStyle',':'); 
plot(freqV(plS),PV(plS),'r','LineStyle',':');
hold on
plot(freqE,PE,'LineStyle',':'); 
plot(freqE(plS),PE(plS),'r','LineStyle',':');

ES = sum(PS(plS))/sum(PS)
EA = sum(PA(plS))/sum(PA)
EV = sum(PV(plS))/sum(PV)
EE = sum(PE(plS))/sum(PE)