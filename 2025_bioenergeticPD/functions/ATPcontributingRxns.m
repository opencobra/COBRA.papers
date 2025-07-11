
dir = 1; % 0, consumption rxns
transportRxns = {'ATPtm','ATPtn','ATPtx','ATP1ter','ATP2ter','EX_atp[e]','DNDPt13m','DNDPt2m','DNDPt31m'...
    ,'DNDPt56m','DNDPt32m','DNDPt57m','DNDPt20m','DNDPt44m','DNDPt19m','DNDPt43m','r1116'};
met2test = {'atp[c]','atp[m]','atp[n]','atp[r]','atp[x]', ...
    'atp[l]', 'atp[i]', 'atp[e]', 'atp[p]'}; % find atp related rxns
tmpModel.mets = model.mets;
isIncluded = ~ismember(model.rxns,transportRxns);
tmpModel.S = model.S(:,isIncluded);
tmpV = FBAsolution_maxATP1.v(isIncluded);
[P,C,vP,vC] = computeFluxSplits(tmpModel,met2test,tmpV,1); %using only the sign of the stoichiometric coefficient (coeffSign =1).

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

% collect results for model
metprod = find(vMetAll);
RxnsNamesAll = model.rxns(metprod);
vMetAll = vMetAll(metprod);
metprod_phi = metprod_phi(metprod);
metRs = [RxnsNamesAll num2cell([vMetAll metprod_phi])];