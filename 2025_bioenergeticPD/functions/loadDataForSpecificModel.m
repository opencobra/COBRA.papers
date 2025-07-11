function [coreRxnAbbr,coreRxnWeights,noRxnAbbr,noRxnWeights,EntrezGeneID,EntrezGeneIDWeights] = loadDataForSpecificModel(modelGeneric,fileName)
% read in data for creation of a cell/tissue/organ specific model
%
%INPUT
% modelGeneirc      Model in COBRA toolbox format e.g. ReconX
% filename   xls spreadsheet with the following tabs
%             activeReactions     reactions that must be part of the core set
%             inactiveReactions   reactions that must not be part of the core set
%             genes               genes that are expressed within the cell (optionally with weights)
%             transportGenes      genes for transporters that are expressed within the cell
%
%OPTIONAL INPUT
% fluxEpsilon   {1e-4} lowest value of flux that is considered non-zero in
%               flux consistency check
%
%
%OUTPUT
% coreRxnAbbr   reactions that must be in the specific model
% noRxnAbbr     reactions that must not be in the specific model
%
% Ines Thiele, June 2015
%
% Longfei Mao, 4/11/2015   filanme can be a structure variable produced by 
%                          the "postExpression" function

%Read in the data from each tab of the input file
%may need to do this with tab delimited files instead of xls to ensure
%cross platform compatibility

%TODO - setup input files and createSpecificModel.m to be able to handle
%weights on reactions (and genes)
coreRxnWeights=[];
noRxnWeights=[];

if isstruct(fileName) % check if the fileName is a variable structure
    omicsDataStru=fileName;
    if ~isempty(omicsDataStru.activeReactions)
        RxnAbbr=omicsDataStru.activeReactions(:,1); % 'activeReactions'
    else 
        RxnAbbr= {};
    end
        coreRxnAbbr = RxnAbbr;
    try
        noRxnAbbr = omicsDataStru.inactiveReactions; % 'inactiveReactions
    catch
        disp('good');
    end
else

    %load reactions known to be carried out by the cell/tissue/organ
    %these reactions shall be added to the core reaction set
    
    [Numbers, Strings] = xlsread(fileName,'activeReactions');
    
    if (~isempty(Strings))
        RxnAbbr = Strings(2:end,2);
    end
    
    %load reactions known NOT to be carried out by the cell/tissue/organ
    %these reactions should NOT be in the cell/tissue/organ specific model
    [Numbers, Strings] = xlsread(fileName,'inactiveReactions');
    if (~isempty(Strings))
        noRxnAbbr = Strings(2:end,2);
    end
    
    % load in gene weighting data for cell/tissue/organ
end

if isstruct(fileName)
    Numbers=str2double(omicsDataStru.genes(:,1));          % 'genes'
    %coreRxnWeights=omicsDataStru.genes(:,2)
else
    [Numbers, Strings] = xlsread(fileName,'genes');
end


if (~isempty(Numbers))
    EntrezGeneID = Numbers(:,1);
    % add reactions from HPA to core set
    EntrezGeneID(isnan(EntrezGeneID))=[];
    EntrezGeneID = cellstr(num2str(EntrezGeneID));
    for i = 1 : length(EntrezGeneID)
        EntrezGeneID(i) = regexprep(EntrezGeneID(i),' ','');
    end
    %TODO - set up script to read the weights
    EntrezGeneIDWeights=ones(length(EntrezGeneID),1);
else
    EntrezGeneID = '';
end



% if isstruct(fileName)
%     Numbers=omicsDataStru.transportGenes;
%     Numbers=str2double(Numbers);
% else
%     %transporter genes that have been reported to be expressed (and active) in the cell/tissue/organ
%     %all reactions associated with these genes will be in the core set
%     [Numbers, Strings] = xlsread(fileName,'transportGenes');
% end
% 
% if (~isempty(Numbers))
%     Trgenes = Numbers(:,1);
%     % add reactions from transport to core set
%     Trgenes(isnan(Trgenes))=[];
%     Trgenes = cellstr(num2str(Trgenes));
%     for i = 1 : length(Trgenes)
%         Trgenes(i) = regexprep(Trgenes(i),' ','');
%     end
%     [TrReactionList] = getRxnsFromGenes(modelGeneric, Trgenes,0);
%     
%     if 0
%         %TODO - this can be vectorised
%         
%         % add all exchange reactions for those transport metabolites to Core list as well
%         sa = printRxnFormula(model,TrReactionList);
%         
%         %list of exhange reactions
%         modelexchanges1 = strmatch('Ex_',model.rxns);
%         modelexchanges2 = strmatch('EX_',model.rxns);
%         modelexchanges = [model.rxns(modelexchanges1);model.rxns(modelexchanges2)];
%         CoreEx =[];
%         for i = 1 : length(sa)
%             [metaboliteList,stoichCoeffList,revFlag] = parseRxnFormula(sa{i});
%             % grab all reactions that this metabolite appears in
%             [rxnList, rxnFormulaList] = findRxnsFromMets(model, metaboliteList);
%             % grab corresponding exchange reaction
%             CoreEx = [CoreEx;intersect(modelexchanges,rxnList)];
%         end
%         CoreEx = unique(CoreEx);
%         coreRxns = [RxnAbbr;TrReactionList;CoreEx];
%     else
%         coreRxnAbbr = [RxnAbbr;TrReactionList];
%     end
% else
%     coreRxnAbbr = RxnAbbr;
% end

coreRxnAbbr = RxnAbbr;

% add the following additional core reactions by default
%TODO- add these to new table to be read in by this file
cntR = length(coreRxnAbbr)+1;
%coreRxnAbbr(cntR) = {'biomass_maintenance_noTrTr'}; cntR = cntR+1; %'biomass_maintenance_noTrTr' in Recon2.1
%coreRxnAbbr(cntR) = {'biomass_maintenance'}; cntR = cntR+1; %'biomass_maintenance' in Recon2.1
coreRxnAbbr(cntR) = {'EX_h2o[e]'}; cntR = cntR+1;
coreRxnAbbr(cntR) = {'DM_atp_c_'}; cntR = cntR+1;
coreRxnAbbr(cntR) = {'ATPM'}; cntR = cntR+1;
coreRxnAbbr(cntR) = {'EX_o2[e]'}; cntR = cntR+1;
coreRxnAbbr(cntR) = {'EX_co2[e]'}; cntR = cntR+1;
coreRxnAbbr(cntR) = {'EX_glc_D[e]'}; cntR = cntR+1;
coreRxnAbbr(cntR) = {'EX_hco3[e]'}; cntR = cntR+1;

%unique core reactions
coreRxnAbbr = unique(coreRxnAbbr);
%unique no reactions
noRxnAbbr = unique(noRxnAbbr);


%TODO - sanity check to make sure none of the lists of reactions or genes
%conflict with eachother, e.g. rxn both in core and non core list!

