%Compare the effect of predictor choices in estimated strength of response
%to solar cycle forcing
%
% Ted Amdur
% 12/20/19, updated 10/15/20

load predictor_likelihood_10_15.mat;

tsi = categorical(choices(:,1));
ghg = categorical(choices(:,2));
vol = categorical(choices(:,3));
enso = categorical(choices(:,4));

res = bM(:,2);

predC = dataset({res,'res'},{tsi,'tsi'},...
    {ghg,'ghg'},...
    {vol,'vol'},...
    {enso,'enso'});
mdl = fitlm(predC,'res ~ 1 + tsi + ghg + vol + enso');
b = anova(mdl);
predVar = b.SumSq;
varCont = predVar./(sum(predVar));


%mdl = fitlm(Model_Year,MPG,'CategoricalVars',1,'VarNames',{'Model_Year','MPG'})