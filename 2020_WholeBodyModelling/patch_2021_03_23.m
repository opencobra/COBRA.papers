% patch for Harvey and Harvetta version 1.03c resulting into version 1.03d

male = loadPSCMfile('Harvey');
female = loadPSCMfile('Harvetta');

%% Fix the bounds on sIEC_EX_met(e)_[bc]
% Male
EX = male.rxns(contains(male.rxns,'sIEC_EX_'));
E = male.rxns(contains(male.rxns,'(e)'));
EXE = intersect(EX,E);
% the reactions are written as 
% sIEC_EX_34dhphe(e)_[bc]	sIEC_34dhphe[e] 	->	34dhphe[bc] 
% hence the upper bounds on these reactions should be 0.
male.ub(ismember(male.rxns,EXE)) = 0;

% Fenale
EX = female.rxns(contains(female.rxns,'sIEC_EX_'));
E = female.rxns(contains(female.rxns,'(e)'));
EXE = intersect(EX,E);
% the reactions are written as 
% sIEC_EX_34dhphe(e)_[bc]	sIEC_34dhphe[e] 	->	34dhphe[bc] 
% hence the upper bounds on these reactions should be 0.
female.ub(ismember(female.rxns,EXE)) = 0;

male.version = 'Harvey_1.03d, generated with patch_2021_03_23';
female.version = 'Harvetta_1.03d, generated with patch_2021_03_23';

save Harvey_1_03d male
save Harvetta_1_03d female