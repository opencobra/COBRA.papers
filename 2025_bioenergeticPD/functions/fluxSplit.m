% find [m] rxns
Mrxns=findRxnFromCompartment(ComplexImodel,'[m]');
for n=1:length(Mrxns)
% Mflux{n}=FBA_ComplexI.v(ismember(ComplexImodel.rxns,Mrxns{n}));
Mflux{n}=solutionConsistency.vThermo(ismember(ComplexImodel.rxns,Mrxns{n}));
end
Mflux=Mflux';
submodel=extractSubNetwork(ComplexImodel,Mrxns);

Mflux=cell2mat(Mflux);

submodel=GeneratedModel;
allflux =FBA_ComplexI.v;

met2test = {'h[c]'};

% ATPprod = {'ATPS4mi','PGK','PYK','SUCOASm'}';

transportRxns = {'ATPtm','DNDPt13m'...
    ,'DNDPt2m','DNDPt31m','DNDPt56m','DNDPt32m','DNDPt57m','DNDPt20m','DNDPt44m'...
    ,'DNDPt19m','DNDPt43m','r1116'};
  %% Compute flux splits
    % Remove excluded reactions (transportRxns)
    tmpModel.mets = submodel.mets;
    % isIncluded = ~ismember(submodel.rxns,transportRxns);
    % tmpModel.S = submodel.S(:,isIncluded);
    tmpModel.S=submodel.S;
    % tmpV = Mflux(isIncluded);
    tmpV = allflux;
    [P,C,vP,vC] = computeFluxSplits(tmpModel,met2test,tmpV,1); 
     % decide if production (1) or consumption ~1.
    vMetAll = zeros(size(isIncluded));
    metprod_phi = zeros(size(isIncluded));
    if dir == 1
        vMetAll(isIncluded) = vP;
        metprod_phi(isIncluded) = P;
    else
        vMetAll(isIncluded) = vC;
        metprod_phi(isIncluded) = C;
    end
    Mflux(abs(Mflux)<0.001)=0;
    fluxValue=[Mrxns(:,1),table(Mflux)];
    fluxValue.Properties.VariableNames={'reaction','value'};
    writetable(fluxValue,'~/Documents/rxns_c1.csv')
%%%%%%%%%%%%%%
    tmpModel.mets = ComplexImodel.mets;
    isIncluded = ~ismember(ComplexImodel.rxns,transportRxns);
    tmpModel.S = ComplexImodel.S(:,isIncluded);
    tmpV = solutionConsistency.vThermo(isIncluded);
    [P,C,vP,vC] = computeFluxSplits(tmpModel,met2test,tmpV,1); 
     % decide if production (1) or consumption ~1.
    vMetAll = zeros(size(isIncluded));
    metprod_phi = zeros(size(isIncluded));
    if dir == 1
        vMetAll(isIncluded) = vP;
        metprod_phi(isIncluded) = P;
    else
        vMetAll(isIncluded) = vC;
        metprod_phi(isIncluded) = C;
    end
    %%
   [minFlux, maxFlux, ~, ~] = fluxVariability(ComplexImodel);
   osenseStr='max';
   printLevel=1;
   [vSparse, sparseRxnBool, essentialRxnBool]  = sparseFBA(ComplexImodel, osenseStr, 1, 1, 'cappedL1', printLevel) 
    %Replace the lower and upper bounds by the min and max Fluxes obtained
    %through FVA:
    FBA_ComplexI.v=minFlux;
    FBA_ComplexI.v=maxFlux;
       
    ComplexImodel=model2escher(ComplexImodel);    
    flux=FBA_ComplexI.v;
    fluxValue=[ComplexImodel.rxns,table(flux)];
    fluxValue.Properties.VariableNames={'reaction','value'};
    writetable(fluxValue,'~/Documents/test_TCA.csv')