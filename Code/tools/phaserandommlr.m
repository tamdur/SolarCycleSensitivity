function range = phaserandommlr(series, pred, interval,iter)
%PHASERANDOMMLR Produce lower and upper bounds of confidence interval for MLR
%using iter instances of phase randomization
%   INPUTS:
%       series: NX1 temperature series
%       pred: NXM array of M predictors of observations N, where first
%           column is ones and second column is TSI. Lags should already be
%           included
%       interval: Size of CI. Enter 0.95 for 95%
%       iter: number of iterations of phase randomization
%   OUTPUTS:
%       range: 2X1 vector of lower bound and upper bound for CI

bR = NaN(iter,1);
tsiSub = pred(:,2);
tsiSub(isnan(tsiSub)) = nanmean(pred(:,2));
for ii = 1:iter
    pred(:,2) = phaserandomize(tsiSub);
    b = regress(series,pred);
    bR(ii) = b(2);
end
range(1) = quantile(bR,0.5-interval/2);
range(2) = quantile(bR,0.5+interval/2);
end