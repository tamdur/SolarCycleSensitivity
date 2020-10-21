%Use calcmlrsensitivity and calccmdsensitivity to generate 46-yr
%sensitivity estimates used in Fig. 4 of ASH20
%
% Ted Amdur
% 10/20/20
clearvars
intSt = 1948;
intEnd = 2019;
allSt = 1959; 
allEnd = 2019;
intLen = 46;
stBlock = (intSt:1:(intEnd-intLen+1))'; %2018 is last full year of NCEP record
endBlock = (intSt+intLen-1:1:intEnd)';
%Make CMD calculations
disp('performing rolling CMD calculations');tic;
for ii = 1:length(stBlock)
    ncepCMD(ii) = calccmdsensitivity(stBlock(ii),endBlock(ii),2,1,0,1,2,6/7);
    ncepLDA(ii) = calccmdsensitivity(stBlock(ii),endBlock(ii),2,2,0,1,2,6/7);
    hadCMD(ii) = calccmdsensitivity(stBlock(ii),endBlock(ii),2,1,0,1,1,6/7);
    hadLDA(ii) = calccmdsensitivity(stBlock(ii),endBlock(ii),2,2,0,1,1,6/7);
end
disp('CMD calculations completed');toc;
%Make MLR calculations
disp('performing rolling MLR calculations');tic;
f = loadallforcing(1); %Load available forcing records
for ii = 1:length(stBlock)
    bNCEP = calcmlrsensitivity(stBlock(ii),endBlock(ii),1,1,2,[3 3 2 3 0 0],f);
    bHad = calcmlrsensitivity(stBlock(ii),endBlock(ii),1,1,1,[3 3 2 3 0 0],f);
    ncepMLR(ii) = bNCEP(2);
    hadMLR(ii) = bHad(2);
end
disp('MLR calculations completed');toc;

%Get full interval
[ncepCMDAll,ncepCMDRange] = calccmdsensitivity(allSt,allEnd,2,1,0,1,2,6/7,1,0.95,10000);
[ncepLDAAll, ncepLDARange] = calccmdsensitivity(allSt,allEnd,2,2,0,1,2,6/7,1,0.95,10000);
[hadCMDAll, hadCMDRange] = calccmdsensitivity(allSt,allEnd,2,1,0,1,1,6/7,1,0.95,10000);
[hadLDAAll, hadLDARange] = calccmdsensitivity(allSt,allEnd,2,2,0,1,1,6/7,1,0.95,10000);
[bNCEPAll,ncepRange] = calcmlrsensitivity(allSt,allEnd,1,1,2,[3 3 2 3 0 0],f);
[bHadAll,hadRange] = calcmlrsensitivity(allSt,allEnd,1,1,1,[3 3 2 3 0 0],f);
ncepMLRAll = bNCEPAll(2);
hadMLRAll = bHadAll(2);

save('../Data/code_generated/running_record_ann_21_10.mat','stBlock','endBlock','ncepCMD','ncepLDA',...
        'hadCMD','hadLDA','ncepMLR','hadMLR','ncepCMDAll','ncepCMDRange',...
        'ncepLDAAll','ncepLDARange','hadCMDAll','hadCMDRange','hadLDAAll',...
        'hadLDARange','ncepRange','ncepMLRAll','hadRange','hadMLRAll');


    