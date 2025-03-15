% work on Small inner mitochondrial space model (IMS model)
clear
clc
%%
load('~/fork-cobratoolbox/test/models/mat/Recon3DModel_301.mat');
% model=Recon3D;
OxPRxns={'ATPS4mi','NADH2_u10mi','CYOR_u10mi','CYOOm2i','ATPS4mi'};%,'PIt2mi'
IMSmodel=extractSubNetwork(model,OxPRxns);
IMSmodel.lb(:)=-1000;
IMSmodel.ub(:)= 1000;
% change bound for Oxidative phosphorylation rxns
% IMSmodel.lb(ismember(IMSmodel.rxns,'NADH2_u10mi'))=0.01;
% IMSmodel.lb(ismember(IMSmodel.rxns,'ATPS4mi'))=0.01;
% IMSmodel.lb(ismember(IMSmodel.rxns,'CYOR_u10mi'))=0.01;
% IMSmodel.lb(ismember(IMSmodel.rxns,'CYOOm2i'))=0.01;

% add new rxns for IMS model
dataFolder = '~/drive/bioenergeticsPD/fromXi/reogranizeData/new';
bibliomicData = 'InnerSpaceModelRxns.xlsx'; % revised rxns2constrain file for synaptic_bibliomicData_reconstruction.xlsx
specificData = preprocessingOmicsModel([dataFolder filesep bibliomicData], 1, 1);

% geneRules for rxns2add in correct format
if isfield(specificData, 'rxns2add')
    if ismember('geneRule', specificData.rxns2add.Properties.VariableNames)
        if isnumeric(specificData.rxns2add.geneRule)
            specificData.rxns2add.geneRule = num2cell(specificData.rxns2add.geneRule);
        else
            specificData.rxns2add.geneRule = regexprep(specificData.rxns2add.geneRule, '\.\d', '');
        end
    end
end
param.TolMinBoundary= -1000;
param.TolMaxBoundary= 1000;
param.printLevel=1;
feasTol = getCobraSolverParams('LP', 'feasTol');
%%
% Add reactions (requires: rxns, mass balanced rxnFormulas
% optional:lb, ub, subSystems, grRules to add to the model)
if isfield(specificData, 'rxns2add') && ~isempty(specificData.rxns2add)
    
    if param.printLevel > 0
        disp(' ')
        disp('--------------------------------------------------------------')
        disp(['Adding ' num2str(numel(specificData.rxns2add.rxns)) ' reactions:' ])
        %disp(table([specificData.rxns2add.rxns, specificData.rxns2add.rxnFormulas]))
        disp(' ')
    end
    
    % check if reaction formulas are provided and test if all reactions
    % have assigned a reaction formula
    if ~ismember('rxnFormulas', specificData.rxns2add.Properties.VariableNames)
        disp('No reaction formula provided. Reactions will not be added.');
    elseif any(cellfun('isempty',specificData.rxns2add.rxnFormulas))
        error(['Reaction(s): ' specificData.rxns2add.rxns{cellfun('isempty', specificData.rxns2add.rxnFormulas)} ...
            ' in specificData.rxns2add.rxns requires a reaction formula in specificData.rxns2add.rxnFormulas']);
    else
        
        % Assign subSystem if not provided
        if ~ismember('subSystem', specificData.rxns2add.Properties.VariableNames)
            specificData.rxns2add.subSystem = repelem({'Miscellaneous'}, numel(specificData.rxns2add.rxns), 1);
        end
        % Assign rxnNames (if there is not rxnNames data)
        if ~ismember('rxnNames', specificData.rxns2add.Properties.VariableNames)
            specificData.rxns2add.rxnNames = repelem({'Custom reaction'}, numel(specificData.rxns2add.rxns), 1);
        end
        % Assign grRules (if there is not grRules data)
        if ~ismember('geneRule', specificData.rxns2add.Properties.VariableNames)
            specificData.rxns2add.geneRule = repelem({''}, numel(specificData.rxns2add.rxns), 1);
        end
        % Check if the reactions are already in the model
        if any(ismember(specificData.rxns2add.rxns, IMSmodel.rxns))
            
            if param.printLevel > 0
                disp([num2str(sum(ismember(specificData.rxns2add.rxns, IMSmodel.rxns))) ...
                    ' reaction(s) is(are) already present in the model and will not be added:'])
                disp(specificData.rxns2add(ismember(specificData.rxns2add.rxns, IMSmodel.rxns),:))
            end
            
            % Delete repeated reactions from the input data
            specificData.rxns2add(ismember(specificData.rxns2add.rxns, IMSmodel.rxns),:) = [];
        end
        
        % Add reactions
        for i = 1:length(specificData.rxns2add.rxns)
            IMSmodel = addReaction(IMSmodel, specificData.rxns2add.rxns{i}, 'reactionFormula', ...
                specificData.rxns2add.rxnFormulas{i}, 'subSystem', specificData.rxns2add.subSystem{i},...
                'reactionName', specificData.rxns2add.rxnNames{i}, 'lowerBound', ...
                specificData.rxns2add.lowerBound(i), 'upperBound', specificData.rxns2add.upperBound(i),...
                'geneRule', specificData.rxns2add.geneRule{i}, 'printLevel', param.printLevel);
        end
        %attempts to finds the reactions in the model which export/import from the model
        %boundary i.e. mass unbalanced reactions
        %e.g. Exchange reactions
        %     Demand reactions
        %     Sink reactions
        IMSmodel = findSExRxnInd(IMSmodel, []);
    end    
end
%% check feasibility
rxn1=IMSmodel.rxns(contains(IMSmodel.rxns,'ATPS4mi'));
IMSmodel=changeObjective(IMSmodel,rxn1);
sol = optimizeCbModel(IMSmodel);
if  sol.stat ~= 1
    error('Infeasible stoichiometrically consistent model after adding new reactions.')
else
    if param.printLevel>0
        disp(' ')
        fprintf('%s\n\n','Feasible stoichiometrically consistent model with new reactions.')
    end
end
%% set Thermol model
IMSmodel.metFormulas{14}='C4H4O5';
IMSmodel.metFormulas{15}='C4H2O5'; %Recon3D.metFormulas(ismember(Recon3D.mets,'mal_L[m]'))
% IMSmodel.metFormulas{16}='HO4P';
% IMSmodel.metFormulas{16}='H';
compartments={'m';'i'};
ph=double(string({8,7.2}));
is={'0.15','0.15'};
chi=double(string({'-155','0'}));
T=310.15;
DrGt0_Uncertainty_Cutoff = 20;
confidenceLevel=0.95;
printLevel=2;
concMinDefault=1e-5; % Lower bounds on metabolite concentrations in mol/L
concMaxDefault=1e-2; % Upper bounds on metabolite concentrations in mol/L
model=IMSmodel;

SetThermoModel % use standard met Gibbs energy from Recon3DModel_thermo

IMSmodel=model;

%%
clear param
param.solver = 'mosek';
param.internalNetFluxBounds='original';
param.externalNetFluxBounds='original';
param.method='fluxes';
param.printLevel=2;
param.debug = 1;
param.feasTol=1e-7;
ComplexImodel=IMSmodel;
% rxn2=ComplexImodel.rxns(contains(ComplexImodel.rxns,'ATPS4mi'));
rxn2=ComplexImodel.rxns(contains(ComplexImodel.rxns,'DM_atp[m]'));
ComplexImodel = changeObjective(ComplexImodel,rxn2);

if ~isfield(ComplexImodel,'SConsistentRxnBool')
massBalanceCheck=1;
[SConsistentMetBool, SConsistentRxnBool, SInConsistentMetBool, SInConsistentRxnBool, unknownSConsistencyMetBool, unknownSConsistencyRxnBool, ComplexImodel, stoichConsistModel] = findStoichConsistentSubset(ComplexImodel, ...
    massBalanceCheck);
end
ComplexImodel=stoichConsistModel;

% change cf/cr
N = ComplexImodel.S(:,ComplexImodel.SConsistentRxnBool);
[m,n] = size(N);
ComplexImodel.cf=zeros(n,1);
ComplexImodel.cr=zeros(n,1);
ComplexImodel.cf = ComplexImodel.DrGt0(ComplexImodel.SConsistentRxnBool);
ComplexImodel.cr=ComplexImodel.cf;
% ComplexImodel.cr = abs(ComplexImodel.DrGt0(ComplexImodel.SConsistentRxnBool));
% % ComplexImodel.cf = -abs(N'* IMSmodel.DfGt0);
% ComplexImodel.cf=-ComplexImodel.cr;
% ind=find(ismember(ComplexImodel.rxns(ComplexImodel.SConsistentRxnBool),'Htmi'));
ind=find(ismember(ComplexImodel.rxns(ComplexImodel.SConsistentRxnBool),'MDHm'));
ComplexImodel.cf(ind)=-1;
ComplexImodel.cr(ind)=6;

FBA_ComplexI=entropicFluxBalanceAnalysis(ComplexImodel,param);

%%
% chemical potential y_N
N=ComplexImodel.S(find(ComplexImodel.SIntMetBool),find(ComplexImodel.SConsistentRxnBool));
rxns=ComplexImodel.rxns(ComplexImodel.SConsistentRxnBool);
potential=N'*FBA_ComplexI.y_N; % chemical potential for all the internal rxns
InternalPotential=[table(rxns),table( potential)];

% for solution part
% [PR, PN, PC, PL] = subspaceProjector(N, 0, 'all');
% SolutionPotential=PR'*FBA_ComplexI.y_N; % chemical potential for solution 

% % y_N= yL + yC; N'*y_N = N'*yL + N'*yC;  N'*yL = 0;
% [yR, yN, yC, yL] = subspaceProjector(FBA_ComplexI.y_N, 0, 'all');
% Null=N'*yL; % not all zero??
% SolutionPotential=N'*yC;

% save flux as rxns data, chemical potential as mets data for Escher map
fluxValue=[ComplexImodel.rxns,table(FBA_ComplexI.v)];
fluxValue.Properties.VariableNames={'reaction','value'};
rxnpotential=[ComplexImodel.rxns(ComplexImodel.SIntRxnBool),table(potential)];
rxnpotential.Properties.VariableNames={'gene','value'};% use rxnspotential as gene data
metpotental=[ComplexImodel.mets,table(FBA_ComplexI.y_N)];
metpotental.Properties.VariableNames={'metabolites','value'};
fluxValue.value(abs(fluxValue.value)<0.0001)=0;
rxnpotential.value(abs(rxnpotential.value)<0.0001)=0;
metpotental.value(abs(metpotental.value)<0.0001)=0;
% writetable(fluxValue,'~/Documents/innerSpaceFlux.csv')
% writetable(rxnpotential,'~/Documents/innerSpaceRxnsPotential.csv')
% writetable(metpotental,'~/Documents/innerSpaceMetsPotential.csv')
potential(ismember(rxns,'PIt2mi'))
potential(ismember(rxns,'NADH2_u10mi'))
potential(ismember(rxns,'SUCD1m'))
potential(ismember(rxns,'CYOR_u10mi'))
potential(ismember(rxns,'CYOOm2i'))
potential(ismember(rxns,'ATPS4mi'))
potential(ismember(rxns,'Htmi'))

% the chemical potential of h[i] could be higher if we force to consume
% nadh[m]->(-1000,-1) or generate nad[m]->(1,1000)
%% check IMSmodel: log(Vr/Vf) = N'*y_N + (Zvf - Zvr) + Z_vi = model.DrGt0?
ind=find(ismember(ComplexImodel.rxns(ComplexImodel.SIntRxnBool&ComplexImodel.SConsistentRxnBool),{'ATPS4mi'}));%'Htmi','ATPS4mi'

logVrVf=log(FBA_ComplexI.vr./FBA_ComplexI.vf);
logVrVf(ind)

N=ComplexImodel.S(find(ComplexImodel.SIntMetBool),find(ComplexImodel.SIntRxnBool& ComplexImodel.SConsistentRxnBool));
potential=N'*FBA_ComplexI.y_N;
potential(ind)

dualZ=FBA_ComplexI.z_vf- FBA_ComplexI.z_vr;
dualZ(ind)

model.DrGt0(ismember(ComplexImodel.rxns,{'Htmi'}))

%% change Reconx to escher format
%IMSmodel1=model2escher(IMSmodel);   
% change to json file using Python
% save '~/Documents/IMSmodel.mat'   IMSmodel -mat
% then run convertToJson.py
%% explore the sum of nad[m] and sum of nadh[m]??
tmpModel.mets = IMSmodel.mets;
isIncluded = ismember(IMSmodel.rxns,IMSmodel.rxns);
tmpModel.S = IMSmodel.S(:,isIncluded);
tmpV = FBA_ComplexI.v(isIncluded);

met2test = {'nad[m]'};
[P1,C1,vP1,vC1] = computeFluxSplits(tmpModel,met2test,tmpV,1); % P= production proportion; vP= production flux; C = consumption proportion; vC= consumption flux
vP1
met2test = {'nadh[m]'};
[P2,C2,vP2,vC2] = computeFluxSplits(tmpModel,met2test,tmpV,1); % P= production proportion; vP= production flux; C = consumption proportion; vC= consumption flux
round(sum(vP2-vC2),4)

%% work on large model
clear
clc
%%
load('~/fork-cobratoolbox/papers/2023_iDopaNeuro/Recon3DModel_301_thermo.mat');
Recon3DThermo=model;
%%
% genericModelName = 'Recon3DModel_301_xomics_input.mat';
% load('~/drive/sbgCloud/projects/xomicsToModel/data/Recon3D_301/Recon3DModel_301_xomics_input.mat');
load('~/drive/bioenergeticsPD/fromXi/model/thermokernel/new/test_SYN_model/test_SYN_model_rawXomicWith[i].mat')
% load('~/drive/bioenergeticsPD/fromXi/model/thermokernel/new/test_model_with[i]/test_model.mat')
model=GeneratedModel;
% model = addReaction(model, 'DM_nad[m]', 'reactionFormula', ...
%             'nad[m] <=>', 'subSystem', 'DM',...
%             'lowerBound', 10, 'upperBound', 1000);
% model = addReaction(model, 'DM_nadh[m]', 'reactionFormula', ...
%             'nadh[m] <=>', 'subSystem', 'DM',...
%             'lowerBound', -1000, 'upperBound', 1000);
% model=removeRxns(model,'NADtm');

% 34HPLFM: atp[m] + h[m] + nmn[m] -> nad[m] + ppi[m]

% model.lb(ismember(model.rxns,'34HPLFM'))=10; % cannot reverse
% model.lb(ismember(model.rxns,'NADtm'))=0;
% model.lb(ismember(model.rxns,'NMNATm'))=2;
% model=removeRxns(model,'NADtm');
% model = addReaction(model, 'NADHtm', 'reactionFormula', ...
%             'nadh[m] <=> nadh[c]', 'subSystem', 'Transport, mitochondrial',...
%             'lowerBound', -1000, 'upperBound', 1000); 
%% set Thermol model
% for SYN model
model.metFormulas(ismember(model.mets,'CE1554[c]'))={'C5H8NO3'};
model.metFormulas(ismember(model.mets,'CE1554[m]'))={'C5H8NO3'};
model.metCompartments(ismember(model.mets,'CE1554[c]'))={'c'};
model.metCompartments(ismember(model.mets,'CE1554[m]'))={'m'};
model=rmfield(model,{'concMin','concMax','DfGt0','DrGt0'});
% use standard met Gibbs energy from Recon3DThermo, new mets have zero 
model.metFormulas(ismember(model.mets,'etfqh2[m]'))={'R'};
model.metFormulas(ismember(model.mets,'etfq[m]'))={'R'};

T=model.T;
compartments=model.compartments(1:9);
ph=model.ph(1:9);
is=model.is(1:9);
chi=model.chi(1:9);
concMinDefault=1e-5; % Lower bounds on metabolite concentrations in mol/L
concMaxDefault=1e-2; % Upper bounds on metabolite concentrations in mol/L
confidenceLevel=0.95;
%%
SetThermoModel
%%
clear param
param.solver = 'mosek';
param.internalNetFluxBounds='original';
param.externalNetFluxBounds='original';
param.method='fluxes';
param.printLevel=2;
param.debug = 1;
param.feasTol=1e-7;
ComplexImodel=model;
% rxn2=ComplexImodel.rxns(contains(ComplexImodel.rxns,'ATPS4mi'));
if ~ismember('ATPM',ComplexImodel.rxns) &  ismember('DM_atp_c_',ComplexImodel.rxns)
   ComplexImodel.rxns(ismember(ComplexImodel.rxns,'DM_atp_c_'))={'ATPM'};
end
rxn2=ComplexImodel.rxns(contains(ComplexImodel.rxns,'ATPM'));

ComplexImodel = changeObjective(ComplexImodel,rxn2);

% %change cf/cr
% N = ComplexImodel.S(:,ComplexImodel.SConsistentRxnBool);
% [m,n] = size(N);
% ComplexImodel.cf=zeros(n,1);
% ComplexImodel.cr=zeros(n,1);
% % ComplexImodel.cr = abs(ComplexImodel.DrGt0(ComplexImodel.SConsistentRxnBool));
% ComplexImodel.cf=-abs(ComplexImodel.DrGt0(ComplexImodel.SConsistentRxnBool));
% ComplexImodel.cr=abs(ComplexImodel.DrGt0(ComplexImodel.SConsistentRxnBool));

% ComplexImodel.cf(ismember(ComplexImodel.rxns(ComplexImodel.SConsistentRxnBool),'Htmi'))=0.1;
% ComplexImodel.cr(ismember(ComplexImodel.rxns(ComplexImodel.SConsistentRxnBool),'Htmi'))=15;


FBA_ComplexI=entropicFluxBalanceAnalysis(ComplexImodel,param);
FBA_ComplexI.v(abs(FBA_ComplexI.v)<1e-4)=0;

%% check log(Vr/Vf) = N'*y_N + (Zvf - Zvr) + Z_vi = model.DrGt0?
ind=find(ismember(ComplexImodel.rxns(ComplexImodel.SIntRxnBool&ComplexImodel.SConsistentRxnBool),{'PIt2mi'}));%'Htmi','ATPS4mi'

logVrVf=log(FBA_ComplexI.vr./FBA_ComplexI.vf);
logVrVf(ind)

N=ComplexImodel.S(find(ComplexImodel.SIntMetBool),find(ComplexImodel.SIntRxnBool& ComplexImodel.SConsistentRxnBool));
potential=N'*FBA_ComplexI.y_N;
potential(ind)

FBA_ComplexI.z_vi(ind)

dualZ=FBA_ComplexI.z_vf- FBA_ComplexI.z_vr;
dualZ(ind)

FBA_ComplexI.v(ismember(ComplexImodel.rxns,{'10FTHFtm'}))

model.DrGt0(ismember(ComplexImodel.rxns,{'PPItm'}))

%
rxns={{'NADH2_u10mi'},}
x=sum(ComplexImodel.SIntRxnBool&ComplexImodel.SConsistentRxnBool);
ComplexImodel=findSExRxnInd(ComplexImodel)
y=FBA_ComplexI.v(ComplexImodel.SIntRxnBool&ComplexImodel.SConsistentRxnBool);
plot(1:x,y)
hold on
y1=FBA_ComplexI.vf;
scatter(1:x,y1)
y2=FBA_ComplexI.vr;
scatter(1:x,y2)
y3=FBA_ComplexI.vf-FBA_ComplexI.vr;
scatter(1:x,y3)
legend('v','vf','vr','vf-vr')
InternalRxns=ComplexImodel.rxns(ComplexImodel.SIntRxnBool&ComplexImodel.SConsistentRxnBool);
All=table(InternalRxns,y,y1, y2,y3)
All.Properties.VariableNames={'rxns','v','vf','vr','vf-vr'}
%%        
N=ComplexImodel.S(find(ComplexImodel.SIntMetBool),find(ComplexImodel.SIntRxnBool));
rxns=ComplexImodel.rxns(ComplexImodel.SIntRxnBool);
[yR, yN, yC, yL] = subspaceProjector(FBA_ComplexI.y_N, 0, 'all');
% potential=N'*yC; 
potential=N'*FBA_ComplexI.y_N;
fluxValue=[ComplexImodel.rxns,table(FBA_ComplexI.v)];
fluxValue.Properties.VariableNames={'reaction','value'};
% writetable(fluxValue,'~/Documents/test_TCA.csv')


rxnpotential=[ComplexImodel.rxns(ComplexImodel.SIntRxnBool),table(potential)];
rxnpotential.Properties.VariableNames={'gene','value'};% use rxnspotential as gene data
mets=table(ComplexImodel.mets);
mets.Properties.VariableNames={'metabolites'}
metpotental=[mets,table(FBA_ComplexI.y_N)];
metpotental.Properties.VariableNames={'metabolites','value'};
fluxValue.value(abs(fluxValue.value)<0.0001)=0;
rxnpotential.value(abs(rxnpotential.value)<0.0001)=0;
metpotental.value=round(metpotental.value,4);
metpotental.value(abs(metpotental.value)<0.0001)=0;

rxnpotential(ismember(model.rxns(model.SIntRxnBool),'Htmi'),:) % should be a negative potential
metpotental(ismember(model.mets,{'h[i]'}),:)
metpotental(ismember(model.mets,{'h[m]'}),:)
% rxnpotential(ismember(model.rxns(model.SIntRxnBool),rxns),:)
% fluxValue(ismember(model.rxns,rxns),:)
% printRxnFormula(model, rxns);

% mets={'h[i]','h[m]','nad[m]','nadh[m]','q10[m]','q10h2[m]'};
% metpotental(ismember(model.mets(model.SIntMetBool),mets),:)
% 
% rxns=model.rxns(find(model.S(ismember(model.mets,'nad[m]'),:)>0));
% rxnsfomula=printRxnFormula(model, rxns);
% NADgeneration=[string(rxns) rxnsfomula model.lb(ismember(model.rxns,rxns)) model.ub(ismember(model.rxns,rxns)),model.subSystems(ismember(model.rxns,rxns)) FBA_ComplexI.v(ismember(model.rxns,rxns))];
% rxns=model.rxns(find(model.S(ismember(model.mets,'nad[m]'),:)<0));
% rxnsfomula=printRxnFormula(model, rxns);
% NADconsumption=[string(rxns) rxnsfomula model.lb(ismember(model.rxns,rxns)) model.ub(ismember(model.rxns,rxns)),model.subSystems(ismember(model.rxns,rxns)) FBA_ComplexI.v(ismember(model.rxns,rxns))];
% 
% NADgeneration=sortrows(NADgeneration,5);
% NADconsumption=sortrows(NADconsumption,5);
% 
% model.lb(ismember(model.rxns,NADconsumption(:,1)))=0;
%
%% explore the sum of nad[m] and sum of nadh[m]??
tmpModel.mets = model.mets;
isIncluded = ismember(model.rxns,model.rxns);
tmpModel.S = model.S(:,isIncluded);
tmpV = FBA_ComplexI.v(isIncluded);

met2test = {'nad[m]'};
[P1,C1,vP1,vC1] = computeFluxSplits(tmpModel,met2test,tmpV,1); % P= production proportion; vP= production flux; C = consumption proportion; vC= consumption flux
sum(vC1)
met2test = {'nadh[m]'};
[P2,C2,vP2,vC2] = computeFluxSplits(tmpModel,met2test,tmpV,1); % P= production proportion; vP= production flux; C = consumption proportion; vC= consumption flux
sum(vP2)
round(sum(vP2-vC1),4)
% change nad/nadh bound
nadhPrxns=model.rxns(vP2~=0); 
% model.lb(ismember(model.rxns,nadhPrxns))=1;
nadCrxns=model.rxns(vC1~=0); 
% model.ub(ismember(model.rxns,nadCrxns))=-0.001;