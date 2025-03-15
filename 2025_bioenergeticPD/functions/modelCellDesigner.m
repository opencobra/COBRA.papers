% Compare the SYN model (in vivo) and iDopaNeuroCT model (in vitro)

% load SYN model (in vivo)
load('~/drive/bioenergeticsPD/fromXi/model/thermokernel/new/test_SYN_model/test_SYN_model.mat')
SYN=GeneratedModel;
% load iDopaNeuroCT/C model (in vitro)
load('~/drive/sbgCloud/data/models/published/iDopaNeuro/iDopaNeuroCT.mat')
load('~/drive/sbgCloud/data/models/published/iDopaNeuro/iDopaNeuroC.mat')

multimodels.SYN=SYN;
multimodels.iDopaNeuroC=iDopaNeuroC;
multimodels.iDopaNeuroCT=iDopaNeuroCT;
%% model comparison
[overlapResults,statistic]=compareXomicsModels(multimodels);
plotOverlapResults(overlapResults,statistic);

%% load iDopaNeuonCombo xmlmap (cellDesigner map)
[xmlModel, mapModel] = transformXML2Map('~/drive/metaPD/data/DN_models/xmlMap/iNESC2DN_ComboMap_7Mar2024');
%% change map with different color or size
mapModelchange = unifyMetabolicMapCD(mapModel);
% not all the overlapped mets in the overlapped rxns.
overlapmets=mapModelchange.specName(ismember(mapModelchange.specName,SYN.mets));
overlaprxns=mapModelchange.rxnName(ismember(mapModelchange.rxnName,SYN.rxns));
% mapModelchange = changeNodesArea(mapModelchange, overlapmets, 30, 60);
% mapModelchange = changeMetColor(mapModelchange, overlapmets, 'LIGHTSTEELBLUE');
% mapModelchange=changeRxnColorAndWidth(mapModelchange, mapModelchange.rxnName,'BLACK',2);
% mapModelchange=changeRxnColorAndWidth(mapModelchange, overlaprxns,'BLUE',3);
%% highlight the overlapped rxns/mets
% [annotation_combo] = parseCD('~/drive/metaPD/data/DN_models/xmlMap/iNESC2DN_ComboMap_7Mar2024.xml');
mapModelchange=modifyReactionsMetabolites(mapModelchange, overlaprxns, overlapmets, 'BLUE', 3);
%% save new map
transformMap2XML(xmlModel, mapModelchange, '~/drive/metaPD/data/DN_models/xmlMap/SYNModel_overlap.xml')
%% flux value
clear param
param.solver = 'mosek';
param.internalNetFluxBounds='original';
param.externalNetFluxBounds='original';
param.method='fluxes';
param.printLevel=2;
param.debug = 1;
param.feasTol=1e-7;
ComplexImodel=SYN;
if ~isfield(ComplexImodel,'c') | ~any(ComplexImodel.c)
    rxn2=ComplexImodel.rxns(contains(ComplexImodel.rxns,'ATPS4mi'));
    % rxn2=ComplexImodel.rxns(contains(ComplexImodel.rxns,'ATPM'));
    ComplexImodel = changeObjective(ComplexImodel,rxn2);
end
FBA_ComplexI=entropicFluxBalanceAnalysis(ComplexImodel,param);
%%
[map, flux2, fluxMap] = addFluxFBAdirectionAndColor(mapModelchange, ComplexImodel, FBA_ComplexI);
[newMap] = addNotes(SYN, map);
transformMap2XML(xmlModel, newMap, '~/drive/metaPD/data/DN_models/xmlMap/SYNModel_eFBAflux.xml')
