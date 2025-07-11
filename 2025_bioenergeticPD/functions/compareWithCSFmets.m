%
clear param
param.solver = 'mosek';
param.internalNetFluxBounds='original';
param.externalNetFluxBounds='original';
param.method='fluxes';
param.printLevel=2;
param.debug = 1;
param.feasTol=1e-7;
% SYN
ComplexImodel=Allmodels.SYN.SYN;
if ~isfield(ComplexImodel,'c') | ~any(ComplexImodel.c)
%     rxn2=ComplexImodel.rxns(contains(ComplexImodel.rxns,'ATPS4mi'));
    rxn2=ComplexImodel.rxns(contains(ComplexImodel.rxns,'ATPM'));
    ComplexImodel = changeObjective(ComplexImodel,rxn2);
%     ind=find(contains(ComplexImodel.rxns,'ICDHxm'));
%     ComplexImodel.c(ind)=1;
end
FBA_ComplexI_SYN=entropicFluxBalanceAnalysis(ComplexImodel,param);
% SYNPD
ComplexImodel=Allmodels.SYNPD.SYNPD;
if ~isfield(ComplexImodel,'c') | ~any(ComplexImodel.c)
%     rxn2=ComplexImodel.rxns(contains(ComplexImodel.rxns,'ATPS4mi'));
    rxn2=ComplexImodel.rxns(contains(ComplexImodel.rxns,'ATPM'));
    ComplexImodel = changeObjective(ComplexImodel,rxn2);
%     ind=find(contains(ComplexImodel.rxns,'ICDHxm'));
%     ComplexImodel.c(ind)=1;
end
FBA_ComplexI_SYNPD=entropicFluxBalanceAnalysis(ComplexImodel,param);
% ASYN
ComplexImodel=Allmodels.ASYN.ASYN;
if ~isfield(ComplexImodel,'c') | ~any(ComplexImodel.c)
%     rxn2=ComplexImodel.rxns(contains(ComplexImodel.rxns,'ATPS4mi'));
    rxn2=ComplexImodel.rxns(contains(ComplexImodel.rxns,'ATPM'));
    ComplexImodel = changeObjective(ComplexImodel,rxn2);
%     ind=find(contains(ComplexImodel.rxns,'ICDHxm'));
%     ComplexImodel.c(ind)=1;
end
FBA_ComplexI_ASYN=entropicFluxBalanceAnalysis(ComplexImodel,param);
% ASYNPD
ComplexImodel=Allmodels.ASYNPD.ASYNPD;
if ~isfield(ComplexImodel,'c') | ~any(ComplexImodel.c)
%     rxn2=ComplexImodel.rxns(contains(ComplexImodel.rxns,'ATPS4mi'));
    rxn2=ComplexImodel.rxns(contains(ComplexImodel.rxns,'ATPM'));
    ComplexImodel = changeObjective(ComplexImodel,rxn2);
%     ind=find(contains(ComplexImodel.rxns,'ICDHxm'));
%     ComplexImodel.c(ind)=1;
end
FBA_ComplexI_ASYNPD=entropicFluxBalanceAnalysis(ComplexImodel,param);
% result table
[a,b]=ismember(Allmodels.SYN.SYN.rxns,Rxns);
Rxns=table(Rxns);
Rxns.SYNflux(b(a))=round(FBA_ComplexI_SYN.v(a),4);

[a,b]=ismember(Allmodels.SYNPD.SYNPD.rxns,Rxns.Rxns);
Rxns.SYNPDflux(b(a))=round(FBA_ComplexI_SYNPD.v(a),4);
[a,b]=ismember(Allmodels.ASYN.ASYN.rxns,Rxns.Rxns);
Rxns.ASYNflux(b(a))=round(FBA_ComplexI_ASYN.v(a),4);

[a,b]=ismember(Allmodels.ASYNPD.ASYNPD.rxns,Rxns.Rxns);
Rxns.ASYNPDflux(b(a))=round(FBA_ComplexI_ASYNPD.v(a),4);