% Generates a venn diagram that describes the sources of exchange metabolites 
% in the iDopaNeuro models.
%
% The graph is saved in:
% ~/work/sbgCloud/programReconstruction/projects/exoMetDN/results/codeResults/iDN1/iDopaNeuro1

clear

% define directories
pathSave = ['~' filesep 'work' filesep 'sbgCloud' filesep 'programReconstruction' ...
    filesep 'projects' filesep 'exoMetDN' filesep 'results' filesep 'codeResults' filesep 'iDN1'];

models = {'iDopaNeuroCT'; 'iDopaNeuroC'};
for j = 1:length(models)
    
    % Load model
    load([pathSave filesep models{j} filesep models{j} '.mat'])
    eval(['model = ' models{j} ';'])
    
    % Identify media metabolites, exchange metabolites and metabolites measured
    % with exometabolomic platforms
    mediaRxns = model.XomicsToModelSpecificData.mediaData.rxns;
    exometRxns = model.XomicsToModelSpecificData.exoMet.rxns;
    for i = 1:length(exometRxns)
        if isempty(exometRxns{i})
            exometRxns{i} = num2str(i);
        end
    end
    exchangeRxns = model.rxns(~model.SIntRxnBool);
    
    % Generate diagram
    set7 = numel(intersect(exchangeRxns, intersect(exometRxns, mediaRxns)));
    vennX([numel(setdiff(mediaRxns, [exometRxns; exchangeRxns])) ...
        numel(intersect(mediaRxns, exometRxns))  - set7 ...
        numel(setdiff(exometRxns, [mediaRxns; exchangeRxns]))  ...
        numel(intersect(exchangeRxns, exometRxns))  - set7 ...
        numel(setdiff(exchangeRxns, [exometRxns; mediaRxns]))  ...
        numel(intersect(mediaRxns, exchangeRxns))  - set7 ...
        set7], .01);
    title({'Venn diagram summarising metabolomic', ['measurements (' models{j} ')']}, 'FontSize', 16)
    
    % Save figure
    savefig([pathSave filesep models{j} filesep 'vennMetabolomics_' models{j}])
    saveas(gcf,[pathSave filesep models{j} filesep 'vennMetabolomics_' models{j}],'png')
    saveas(gcf,[pathSave filesep models{j} filesep 'vennMetabolomics_' models{j}],'eps')
    
end
