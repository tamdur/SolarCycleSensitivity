%Use a version of mlrsensitivity with function inputs replacing user inputs
%in order to replicate table 2 outputs
%
% Ted Amdur
% 9/30/20

 %clearvars
 yrEnd = 2019;
stVec = [repmat([1882; 1959],[6 1]); 1959; repmat([1882; 1959],[2 1]); ...
         repmat([1882; 1959], [7 1]); 1959; repmat([1882; 1959],[2 1])];
endVec = ones(36,1).*yrEnd;
ann = [zeros(17,1); ones(19,1)];
detrendTSI = [ones(2,1); zeros(2,1);ones(15,1); zeros(2,1); ones(15,1)];
Tchoice = [ones(12,1);2;3;3;4;4;ones(14,1);2;3;3;4;4];
choiceVec = [3 3 2 3 0 0;
             3 3 2 3 0 0;
             2 1 1 3 0 0;
             2 1 1 3 0 0;
             3 3 2 2 0 0;
             3 3 2 2 0 0;
             3 3 2 3 0 1;
             3 3 2 3 0 1;
             3 3 2 3 1 0;
             3 3 2 3 1 0;
             3 3 2 3 1 1;
             3 3 2 3 1 1;
             3 3 2 3 0 0;
             3 3 2 3 0 0;
             3 3 2 3 0 0;
             3 3 2 3 0 0;
             3 3 2 3 0 0;
             3 3 2 3 0 0;
             3 3 2 3 0 0;
             2 1 1 3 0 0;
             2 1 1 3 0 0;
             3 3 2 2 0 0;
             3 3 2 2 0 0;
             3 3 2 5 0 0;
             3 3 2 5 0 0;
             3 3 2 3 0 1;
             3 3 2 3 0 1;
             3 3 2 3 1 0;
             3 3 2 3 1 0;
             3 3 2 3 1 1;
             3 3 2 3 1 1;
             3 3 2 3 0 0;
             3 3 2 3 0 0;
             3 3 2 3 0 0;
             3 3 2 3 0 0;
             3 3 2 3 0 0];
lagV = NaN(size(choiceVec));
for ii = 18:36
     f = loadallforcing(ann(ii)); %Load available forcing records
    [bT,rnT,lags] = calcmlrsensitivity(stVec(ii),endVec(ii),ann(ii),detrendTSI(ii),Tchoice(ii),choiceVec(ii,:),f);
    b(ii)=bT(2); rn(ii,:)=rnT; lagV(ii,1:length(lags)) = lags;
end
outCell = [{'MONTHLY MLR'};
{['1882-' num2str(yrEnd) ' best: ' num2str(b(1),'%.3f') 'K/(W m^2})' ' +-' num2str(rn(1,2),'%.3f')]};
{['1959-' num2str(yrEnd) ' best: ' num2str(b(2),'%.3f') 'K/(W m^2})' ' +-' num2str(rn(2,2),'%.3f')]};
{['1882-' num2str(yrEnd) ' LR08: ' num2str(b(3),'%.3f') 'K/(W m^2})' ' +-' num2str(rn(3,2),'%.3f')]};
{['1959-' num2str(yrEnd) ' LR08: ' num2str(b(4),'%.3f') 'K/(W m^2})' ' +-' num2str(rn(4,2),'%.3f')]};
{['1882-' num2str(yrEnd) ' 8-year: ' num2str(b(5),'%.3f') 'K/(W m^2})' ' +-' num2str(rn(5,2),'%.3f')]};
{['1959-' num2str(yrEnd) ' 8-year: ' num2str(b(6),'%.3f') 'K/(W m^2})' ' +-' num2str(rn(6,2),'%.3f')]};
{['1882-' num2str(yrEnd) ' IPO: ' num2str(b(7),'%.3f') 'K/(W m^2})' ' +-' num2str(rn(7,2),'%.3f')]};
{['1959-' num2str(yrEnd) ' IPO: ' num2str(b(8),'%.3f') 'K/(W m^2})' ' +-' num2str(rn(8,2),'%.3f')]};
{['1882-' num2str(yrEnd) ' AMO: ' num2str(b(9),'%.3f') 'K/(W m^2})' ' +-' num2str(rn(9,2),'%.3f')]};
{['1959-' num2str(yrEnd) ' AMO: ' num2str(b(10),'%.3f') 'K/(W m^2})' ' +-' num2str(rn(10,2),'%.3f')]};
{['1882-' num2str(yrEnd) ' IPO+AMO: ' num2str(b(11),'%.3f') 'K/(W m^2})' ' +-' num2str(rn(11,2),'%.3f')]};
{['1959-' num2str(yrEnd) ' IPO+AMO: ' num2str(b(12),'%.3f') 'K/(W m^2})' ' +-' num2str(rn(12,2),'%.3f')]};
{['1959-' num2str(yrEnd) ' NCEP: ' num2str(b(13),'%.3f') 'K/(W m^2})' ' +-' num2str(rn(13,2),'%.3f')]};
{['1882-' num2str(yrEnd) ' Berkeley: ' num2str(b(14),'%.3f') 'K/(W m^2})' ' +-' num2str(rn(14,2),'%.3f')]};
{['1959-' num2str(yrEnd) ' Berkeley: ' num2str(b(15),'%.3f') 'K/(W m^2})' ' +-' num2str(rn(15,2),'%.3f')]};
{['1882-' num2str(yrEnd) ' GISS: ' num2str(b(16),'%.3f') 'K/(W m^2})' ' +-' num2str(rn(16,2),'%.3f')]};
{['1959-' num2str(yrEnd) ' GISS: ' num2str(b(17),'%.3f') 'K/(W m^2})' ' +-' num2str(rn(17,2),'%.3f')]};
{'Annual MLR'};
{['1882-' num2str(yrEnd) ' best: ' num2str(b(18),'%.3f') 'K/(W m^2})' ' +-' num2str(rn(18,2),'%.3f')]};
{['1959-' num2str(yrEnd) ' best: ' num2str(b(19),'%.3f') 'K/(W m^2})' ' +-' num2str(rn(19,2),'%.3f')]};
{['1882-' num2str(yrEnd) ' LR08: ' num2str(b(20),'%.3f') 'K/(W m^2})' ' +-' num2str(rn(20,2),'%.3f')]};
{['1959-' num2str(yrEnd) ' LR08: ' num2str(b(21),'%.3f') 'K/(W m^2})' ' +-' num2str(rn(21,2),'%.3f')]};
{['1882-' num2str(yrEnd) ' 8-year: ' num2str(b(22),'%.3f') 'K/(W m^2})' ' +-' num2str(rn(22,2),'%.3f')]};
{['1959-' num2str(yrEnd) ' 8-year: ' num2str(b(23),'%.3f') 'K/(W m^2})' ' +-' num2str(rn(23,2),'%.3f')]};
{['1882-' num2str(yrEnd) ' DJF: ' num2str(b(24),'%.3f') 'K/(W m^2})' ' +-' num2str(rn(24,2),'%.3f')]};
{['1959-' num2str(yrEnd) ' DJF: ' num2str(b(25),'%.3f') 'K/(W m^2})' ' +-' num2str(rn(25,2),'%.3f')]};
{['1882-' num2str(yrEnd) ' IPO: ' num2str(b(26),'%.3f') 'K/(W m^2})' ' +-' num2str(rn(26,2),'%.3f')]};
{['1959-' num2str(yrEnd) ' IPO: ' num2str(b(27),'%.3f') 'K/(W m^2})' ' +-' num2str(rn(27,2),'%.3f')]};
{['1882-' num2str(yrEnd) ' AMO: ' num2str(b(28),'%.3f') 'K/(W m^2})' ' +-' num2str(rn(28,2),'%.3f')]};
{['1959-' num2str(yrEnd) ' AMO: ' num2str(b(29),'%.3f') 'K/(W m^2})' ' +-' num2str(rn(29,2),'%.3f')]};
{['1882-' num2str(yrEnd) ' IPO+AMO: ' num2str(b(30),'%.3f') 'K/(W m^2})' ' +-' num2str(rn(30,2),'%.3f')]};
{['1959-' num2str(yrEnd) ' IPO+AMO: ' num2str(b(31),'%.3f') 'K/(W m^2})' ' +-' num2str(rn(31,2),'%.3f')]};
{['1959-' num2str(yrEnd) ' NCEP: ' num2str(b(32),'%.3f') 'K/(W m^2})' ' +-' num2str(rn(32,2),'%.3f')]};
{['1882-' num2str(yrEnd) ' Berkeley: ' num2str(b(33),'%.3f') 'K/(W m^2})' ' +-' num2str(rn(33,2),'%.3f')]};
{['1959-' num2str(yrEnd) ' Berkeley: ' num2str(b(34),'%.3f') 'K/(W m^2})' ' +-' num2str(rn(34,2),'%.3f')]};
{['1882-' num2str(yrEnd) ' GISS: ' num2str(b(35),'%.3f') 'K/(W m^2})' ' +-' num2str(rn(35,2),'%.3f')]};
{['1959-' num2str(yrEnd) ' GISS: ' num2str(b(36),'%.3f') 'K/(W m^2})' ' +-' num2str(rn(36,2),'%.3f')]}];
%Write results to a text file
% writecell(outTable,'table2outputs.txt','Delimiter','tab') ;
fid = fopen(['../Data/code_generated/table2outputs' num2str(yrEnd) '.txt'],'w');
CT = outCell.';
fprintf(fid,'%s\n', CT{:});
fclose(fid);
