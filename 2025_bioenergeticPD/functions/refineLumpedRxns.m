% 1. genericModel
genericModelName = 'Recon3D_301.mat'; 
switch genericModelName
    case 'Recon3DModel_301.mat'
        load('~/fork-cobratoolbox/test/models/mat/Recon3DModel_301.mat');
        modelOrig=model;
    case 'Recon3D_301.mat'
         load('~/fork-cobratoolbox/test/models/mat/Recon3D_301.mat');
        model=Recon3D;
end

% load unlumped fatty acid rxns
path = '~/drive/metaPD/data/reconstruction_FAD/data';
filename = 'FattyAcidMetabolism.xlsx';
sheetname = 'UnlumpedRxns';
opts = detectImportOptions([path filesep  filename], 'Sheet', sheetname);
opts = setvartype(opts, 'string');
UnlumpedRxns = readtable([path filesep  filename], opts);
for i=1:length(UnlumpedRxns.RxnID)
    rxn=UnlumpedRxns.RxnID(i);
    if ~isempty(rxn) || ~ismissing(rxn)
        if ismember(rxn,model.rxns) % old rxns
            GPR=UnlumpedRxns.GPR(i);
            if ismissing(GPR)
                GPR=model.grRules(ismember(model.rxns,UnlumpedRxns.RxnID{i}));
            end
            model=removeRxns(model,UnlumpedRxns.RxnID{i});
            % change formula
            model = addReaction(model, UnlumpedRxns.RxnID{i}, 'reactionFormula', ...
                UnlumpedRxns.Formula{i}, 'subSystem', UnlumpedRxns.Subsystem{i},...
                'geneRule', GPR{:});
        elseif ~ismember(rxn,model.rxns) % new rxns
            GPR=UnlumpedRxns.GPR(i);
            if ~ismissing(GPR)
                model = addReaction(model, UnlumpedRxns.RxnID{i}, 'reactionFormula', ...
                    UnlumpedRxns.Formula{i}, 'subSystem', UnlumpedRxns.Subsystem{i},...
                    'geneRule', GPR{:});
            else
                model = addReaction(model, UnlumpedRxns.RxnID{i}, 'reactionFormula', ...
                    UnlumpedRxns.Formula{i}, 'subSystem', UnlumpedRxns.Subsystem{i});
            end
        end
    end
end

% remove all the lumped rxxns
filename = 'FattyAcidMetabolism.xlsx';
sheetname = 'lumpedRxns';
opts = detectImportOptions([path filesep  filename], 'Sheet', sheetname);
opts = setvartype(opts, 'string');
lumpedRxns = readtable([path filesep  filename], opts);
removeList=lumpedRxns.lumped_rxn(lumpedRxns.Satus == 'remove');
model = removeRxns(model,char(removeList));

% add new mets formula
filename = 'FattyAcidMetabolism.xlsx';
sheetname = 'newMets';
opts = detectImportOptions([path filesep  filename], 'Sheet', sheetname);
opts = setvartype(opts, 'string');
newMets = readtable([path filesep  filename], opts);
for i=1:length(newMets.metIDs)
    met = newMets.metIDs(i);
    if ismember(met,model.mets)
        [a,b]=ismember(model.mets,met);
        model.metFormulas{a} = newMets.formulas{i};
        model.metInChIString{a} = newMets.inchis(i);
        model.metSmiles{a} = newMets.smiles(i);
        model.metCharges(a) = double(newMets.charges(i));
    end 
end
% save as Recon3D_refined.mat