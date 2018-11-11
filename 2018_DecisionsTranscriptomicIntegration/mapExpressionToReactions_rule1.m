function [expressionRxns, parsedGPR, gene_used] = mapExpressionToReactions_rule1(model, expressionData)                                          
% Determines the expression data associated to each reaction present in
% the model based on the mapping AND=min and OR=max
%
% USAGE:
%    [expressionRxns parsedGPR] = mapExpressionToReactions(model, expressionData) 
%
% INPUTS:
%	model                   model strusture
%	expressionData          mRNA expression data structure
%       .gene               	cell array containing GeneIDs in the same
%                               format as model.genes
%       .value                  Vector containing corresponding expression
%                               value (FPKM/RPKM)
% OUTPUTS:
%   expressionRxns:         reaction expression, corresponding to model.rxns.
%   parsedGPR:              cell matrix containing parsed GPR rule
%
% .. Authors:
%       - Anne Richelle, May 2017 - integration of new extraction methods 

parsedGPR = GPRparser(model);% Extracting GPR data from model

% Find wich genes in expression data are used in the model

[gene_id, gene_expr] = findUsedGenesLevels_modif(model,expressionData);

% Link the gene to the model reactions
[expressionRxns, gene_used] = selectGeneFromGPR_rule1(model, gene_id, gene_expr, parsedGPR);