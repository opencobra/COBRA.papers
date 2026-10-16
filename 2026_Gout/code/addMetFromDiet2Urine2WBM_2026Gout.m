function model = addMetFromDiet2Urine2WBM_2026Gout(model,Metabolite,MetaboliteName,param,reaction,couplingFactor)
% This function adds reaction transport to WBM for a specified metabolite
%
% INPUT
% model             WBM model structure
% Metabolite        Metabolite to be added to WBM
% param             parameters:
%                   param.Diet = 'y'
%                   param.SIabsorp = 'y' % if SI absorption is happening
%                   param.LIabsorp = 'y'
%                   param.BileDuct = 'y'
%                   param.PeripheralOrgan = 'abbr of organ';
%                   param.Sweat = 'y' % requires Skin to be defined in param.PeripheralOrgan
%                   param.BreastMilk = 'y'% requires Breast to be defined in param.PeripheralOrgan
%                   param.StomachUptake = 'y'% requires Stomach to be defined in param.PeripheralOrgan
%                   param.Blood = 'y' % metabolite (e.g., drug) goes straight
%                           into blood stream
%                   param.sink.organs - adds sink reactions to the cytosol
%                           of the listed organs. List, e.g.,
%                           (e.g., param.sinks.organs = {'Heart_';'Muscle_'}; % make
%                           sure that the organ can take up or synthesise the
%                           metabolite (e.g., by defining organ in
%                           param.PeripheralOrgan)
%                   param.sink.biofluids - adds sink reactions to the listed biofluids, e.g.,
%                           (e.g., param.sinks.biofluids = {'[bc]';'[bp]'}; % make
%                           sure that the organ can take up or synthesise the
%                           metabolite (e.g., by defining organ in
%                           param.PeripheralOrgan)
%                   param.Breath = 'y' % = breathing out of volatiles, will go from [LI] to [a]
%                   param.Transdermal - from skin uptake straight into
%                           blood
%                   param.Inhalation - from air via lung into blood
%                   param.MicrobiomeOnly - will not include kidney
%                           reactions
%                   param.Sublingual - goes straight from diet to blood
%                   Buccal works the same way
%
% reaction          Structure containing reaction information per organ.
%                   Note that transport reactions that go into specialised biofluid
%                   compartments (other than from/to [bc]/[e]) must be specified in the
%                   reaction formula or a default reaction is added.
% couplingFactor    Coupling factor of rxns to the corresponding organ biomass
%                           reaction (default: 20000)
%
%
% OUTPUT
% model       expanded WBM model structure
%
% Ines Thiele, Nov 2020 - 2023

if ~exist('coupling','var')
    couplingFactor = 20000; % WBM default
end
% add the reactions for GI tract
space = {' '};
if isfield(param,'Diet') && ~isempty(param.Diet) && contains(lower(param.Diet),'y')
    
    % 1. 'Diet_EX_glc_D[d]'	'Exchange of D-Glucose (Diet)'	'glc_D[d] -->' 'Exchange/demand reaction'
    substr = [Metabolite,'[d]'];
    [model] = addReactionsHH(model, {['Diet_EX_',Metabolite,'[d]']},strcat('Exchange of',space,MetaboliteName,' (Diet)'), {[substr,' -->']}, {''}, {'Exchange/demand reaction'},couplingFactor);
    
    %    if oldStyle % i.e., only GI tract exists
    % 'GI_EX_glc_D[lu]_[d]'	'Exchange of D-Glucose (GI)'	'glc_D[d]  -> glc_D[lu] '   'Transport, biofluid'
    substr = [Metabolite,'[d]'];
    prod = [Metabolite,'[lu]'];
    [model] = addReactionsHH(model, {['GI_EX_',Metabolite,'[lu]_[d]']},strcat('Exchange of',space,MetaboliteName,' (GI)'), {[substr,' --> ',prod]}, {''}, {'Transport, biofluid'},couplingFactor);
    
    % 'SI_EX_glc_D[lu]_[luSI]'	'Exchange of D-Glucose (SI)'	'glc_D[lu]  -> glc_D[luSI] ' 'Transport, biofluid'
    substr = [Metabolite,'[lu]'];
    prod = [Metabolite,'[luSI]'];
    [model] = addReactionsHH(model, {['SI_EX_',Metabolite,'[lu]_[luSI]']},strcat('Exchange of',space,MetaboliteName,' (SI)'), {[substr,' --> ',prod]}, {''}, {'Transport, biofluid'},couplingFactor);
    
    % 'SI_EX_glc_D[luSI]_[luLI]'	'Exchange of D-Glucose (SI)'	'glc_D[luSI]  -> glc_D[luLI] '	'Transport, biofluid'
    substr = [Metabolite,'[luSI]'];
    prod = [Metabolite,'[luLI]'];
    [model] = addReactionsHH(model, {['SI_EX_',Metabolite,'[luSI]_[luLI]']},strcat('Exchange of',space,MetaboliteName,' (SI)'), {[substr,' --> ',prod]}, {''}, {'Transport, biofluid'},couplingFactor);
    
    % 'LI_EX_glc_D[luLI]_[fe]'	'Exchange of D-Glucose (LI)'	'glc_D[luLI]  -> glc_D[fe] '	'Transport, biofluid'
    substr = [Metabolite,'[luLI]'];
    prod = [Metabolite,'[fe]'];
    [model] = addReactionsHH(model, {['LI_EX_',Metabolite,'[luLI]_[fe]']},strcat('Exchange of',space,MetaboliteName,' (LI)'), {[substr,' --> ',prod]}, {''}, {'Transport, biofluid'},couplingFactor);
    
    % 'Excretion_EX_glc_D[fe]'	'Exchange of D-Glucose (Excretion)'	'glc_D[fe]  -> '	'Exchange/demand reaction'
    substr = [Metabolite,'[fe]'];
    [model] = addReactionsHH(model, {['Excretion_EX_',Metabolite,'[fe]']},strcat('Exchange of',space,MetaboliteName,' (Excretion)'), {[substr,' --> ']}, {''}, {'Exchange/demand reaction'},couplingFactor);
    %end
end
%

%
%
%% uptake by small intestine
if  isfield(param,'SIabsorp') && ~isempty(param.SIabsorp) && contains(lower(param.SIabsorp),'y')
    % 'sIEC_EX_glc_D[luI]_[luSI]'	'Exchange of D-Glucose (sIEC)'	'sIEC_glc_D[luI]  <=> glc_D[luSI] '		'Exchange/demand reaction'
   
    substr = ['sIEC_',Metabolite,'[luI]'];
    prod = [Metabolite,'[luSI]'];
    [model] = addReactionsHH(model, {['sIEC_EX_',Metabolite,'[luI]_[luSI]']},strcat('Exchange of',space,MetaboliteName,' (sIEC)'), {[substr,' <==> ',prod]}, {''}, {'Exchange/demand reaction'},couplingFactor);
    
    if exist('reaction','var') && isfield(reaction,'sIEC')
        organ = 'sIEC';
        % find all relevant reactions
        e = reaction.(organ).rxns.abbr(contains(reaction.(organ).rxns.formula,'[luI]'));
        c = reaction.(organ).rxns.abbr(contains(reaction.(organ).rxns.formula,'[c]'));
        rxnT = intersect(e,c);
        if ~isempty(rxnT)
            if (contains(rxnT,'ABCt'))
                % check that indeed an abc efflux
                for r = 1: length(rxnT)
                    x = find(ismember(reaction.(organ).rxns.abbr,rxnT{r}));
                    if ~isempty(contains(reaction.(organ).rxns.formula(x),'atp[c]'))
                        % add to abbr [bc]
                        %     reaction.(organ).rxns.abbr(x)=strcat( reaction.(organ).rxns.abbr(x),'[bc]');
                        reaction.(organ).rxns.abbr(x)=strcat( reaction.(organ).rxns.abbr(x));
                        % default reaction
                        substr = [organ,'_',Metabolite,'[c]'];
                        prod= [organ,'_',Metabolite,'[e]'];
                        Organ_secretion_gene = {''}; % should become an input questions
                        %            [model] = addReactionsHH(model, {[organ,'_',upper(Metabolite),'t[bc]']},strcat('Transport of ',space,MetaboliteName,' from Cytosol to Extracellular Space (blood) (',organ,')'), {[substr,' <==> ',prod]}, {Organ_secretion_gene}, {'Transport, extracellular'},couplingFactor);
                        [model] = addReactionsHH(model, {[organ,'_',upper(Metabolite),'t']},strcat('Transport of ',space,MetaboliteName,' from Cytosol to Extracellular Space (blood) (',organ,')'), {[substr,' <==> ',prod]}, {Organ_secretion_gene}, {'Transport, extracellular'},couplingFactor);
                    end
                end
            end
            model = addOrganRxn(model,reaction, organ,rxnT, couplingFactor,'[luI]');
        else
            % default reaction
            substr = ['sIEC_',Metabolite,'[luI]'];
            prod = ['sIEC_',Metabolite,'[c]'];
            SIEC_uptake_gene = {''}; % should become an input questions
            [model] = addReactionsHH(model, {['sIEC_',upper(Metabolite),'t[luI]']},strcat('Transport of ',space,MetaboliteName,' from Extracellular Space to Cytosol of Mucosal Cells in Small Intestine (sIEC)'), {[substr,' --> ',prod]}, {SIEC_uptake_gene}, {'Transport, extracellular'},couplingFactor);
            
        end
    else
        % default reaction
        substr = ['sIEC_',Metabolite,'[luI]'];
        prod = ['sIEC_',Metabolite,'[c]'];
        SIEC_uptake_gene = {''}; % should become an input questions
        [model] = addReactionsHH(model, {['sIEC_',upper(Metabolite),'t[luI]']},strcat('Transport of ',space,MetaboliteName,' from Extracellular Space to Cytosol of Mucosal Cells in Small Intestine (sIEC)'), {[substr,' --> ',prod]}, {SIEC_uptake_gene}, {'Transport, extracellular'},couplingFactor);
    end
    % transport into the bp
    % e.g., note that gene or mechanism is given here
    % 'sIEC_GLCGLUT2[bpI]'	'transport of glucose into the portal blood (bpI)'	'sIEC_glc_D[c]  -> sIEC_glc_D[bpI] '	'(6514.1)'	'Transport, extracellular'
    if exist('reaction','var') && isfield(reaction,'sIEC')
        organ = 'sIEC';
        % find all relevant reactions
        e = reaction.(organ).rxns.abbr(contains(reaction.(organ).rxns.formula,'[bpI]'));
        c = reaction.(organ).rxns.abbr(contains(reaction.(organ).rxns.formula,'[c]'));
        rxnT = intersect(e,c);
        if ~isempty(rxnT)
            model = addOrganRxn(model,reaction, organ,rxnT, couplingFactor,'[bpI]');
        else
            % default reaction
            substr = ['sIEC_',Metabolite,'[c]'];
            prod = ['sIEC_',Metabolite,'[bpI]'];
            SIEC_secretion_gene = {''}; % should become an input questions
            [model] = addReactionsHH(model, {['sIEC_',upper(Metabolite),'t[bpI]']},strcat('Transport of ',space,MetaboliteName,' into the portal blood (bpI)'), {[substr,' --> ',prod]}, {SIEC_secretion_gene}, {'Transport, extracellular'},couplingFactor);
        end
    else
        % default reaction
        substr = ['sIEC_',Metabolite,'[c]'];
        prod = ['sIEC_',Metabolite,'[bpI]'];
        SIEC_secretion_gene = {''}; % should become an input questions
        [model] = addReactionsHH(model, {['sIEC_',upper(Metabolite),'t[bpI]']},strcat('Transport of ',space,MetaboliteName,' into the portal blood (bpI)'), {[substr,' --> ',prod]}, {SIEC_secretion_gene}, {'Transport, extracellular'},couplingFactor);
    end
    % excretion into bp
    % 'sIEC_EX_glc_D[bpI]_[bp]'	'Exchange of D-Glucose (sIEC)'	'sIEC_glc_D[bpI]  -> glc_D[bp] '	'Transport, biofluid'
    substr = ['sIEC_',Metabolite,'[bpI]'];
    prod = [Metabolite,'[bp]'];
    [model] = addReactionsHH(model, {['sIEC_EX_',Metabolite,'[bpI]_[bp]']},strcat('Transport of ',space,MetaboliteName,' from Extracellular Space to portal vein (sIEC)'), {[substr,' --> ',prod]}, {''}, {'Transport, biofluid'},couplingFactor);
    
end
%% uptake by large intestine
if isfield(param,'LIabsorp') &&~isempty(param.LIabsorp) && contains(lower(param.LIabsorp),'y')
    
    % 'Colon_EX_glc_D[luC]_[luLI]'	'Exchange of D-Glucose (Colon)'	'Colon_glc_D[luC]  <=> glc_D[luLI] '	'Exchange/demand reaction'
    % transport into the colonocytes
    substr = ['Colon_',Metabolite,'[luC]'];
    prod = [Metabolite,'[luLI]'];
    [model] = addReactionsHH(model, {['Colon_EX_',Metabolite,'[luC]_[luLI]']},strcat('Exchange of',space,MetaboliteName,' (Colon)'), {[substr,' <==> ',prod]}, {''}, {'Transport, biofluid'},couplingFactor);
    
    if exist('reaction','var') && isfield(reaction,'Colon')
        organ = 'Colon';
        % find all relevant reactions
        e = reaction.(organ).rxns.abbr(contains(reaction.(organ).rxns.formula,'[luC]'));
        c = reaction.(organ).rxns.abbr(contains(reaction.(organ).rxns.formula,'[c]'));
        rxnT = intersect(e,c);
        if ~isempty(rxnT)
            model = addOrganRxn(model,reaction, organ,rxnT, couplingFactor,'[luC]');
        else
            % default reaction
            % e.g. note that no gene is given here
            % 'Colon_GLCt1r[luC]'	'Glucose Transport (Uniport) (Colon)'	'Colon_glc_D[luC]  <=> Colon_glc_D[c] '		'Transport, extracellular'
            substr = ['Colon_',Metabolite,'[luC]'];
            prod = ['Colon_',Metabolite,'[c]'];
            Colon_uptake_gene = {''}; % should become an input questions
            [model] = addReactionsHH(model, {['Colon_',upper(Metabolite),'t[luC]']},strcat('Transport of ',space,MetaboliteName,' from Extracellular Space to Cytosol of Colonocytes in Large Intestine (Colon)'), {[substr,' <==> ',prod]}, {Colon_uptake_gene}, {'Transport, extracellular'},couplingFactor);
            
        end
    else
        % default reaction
        % e.g. note that no gene is given here
        % 'Colon_GLCt1r[luC]'	'Glucose Transport (Uniport) (Colon)'	'Colon_glc_D[luC]  <=> Colon_glc_D[c] '		'Transport, extracellular'
        substr = ['Colon_',Metabolite,'[luC]'];
        prod = ['Colon_',Metabolite,'[c]'];
        Colon_uptake_gene = {''}; % should become an input questions
        [model] = addReactionsHH(model, {['Colon_',upper(Metabolite),'t[luC]']},strcat('Transport of ',space,MetaboliteName,' from Extracellular Space to Cytosol of Colonocytes in Large Intestine (Colon)'), {[substr,' <==> ',prod]}, {Colon_uptake_gene}, {'Transport, extracellular'},couplingFactor);
    end
    % transport into the bp
    if exist('reaction','var') && isfield(reaction,'Colon')
        organ = 'Colon';
        % find all relevant reactions
        e = reaction.(organ).rxns.abbr(contains(reaction.(organ).rxns.formula,'[bpC]'));
        c = reaction.(organ).rxns.abbr(contains(reaction.(organ).rxns.formula,'[c]'));
        rxnT = intersect(e,c);
        if ~isempty(rxnT)
            model = addOrganRxn(model,reaction, organ,rxnT, couplingFactor,'[bpC]');
        else
            % default reaction
            % e.g., note that no gene is given here
            % 'Colon_GLCt2_2[bpC]'	'Transport of D-Glucose  via Proton Symport (Colon)'	'2 Colon_h[bpC] + Colon_glc_D[bpC]  <=> 2 Colon_h[c] + Colon_glc_D[c] '	'6523.1'	'Transport, extracellular'
            substr = ['Colon_',Metabolite,'[bpC]'];
            prod = ['Colon_',Metabolite,'[c]'];
            Colon_secretion_gene = {''}; % should become an input questions
            [model] = addReactionsHH(model, {['Colon_',upper(Metabolite),'t[bpC]']},strcat('Transport of ',space,MetaboliteName,' into the portal blood (Colon)'), {[substr,' <==> ',prod]}, {Colon_secretion_gene}, {'Transport, extracellular'},couplingFactor);
            
        end
    else
        % default reaction
        % e.g., note that no gene is given here
        % 'Colon_GLCt2_2[bpC]'	'Transport of D-Glucose  via Proton Symport (Colon)'	'2 Colon_h[bpC] + Colon_glc_D[bpC]  <=> 2 Colon_h[c] + Colon_glc_D[c] '	'6523.1'	'Transport, extracellular'
        substr = ['Colon_',Metabolite,'[bpC]'];
        prod = ['Colon_',Metabolite,'[c]'];
        Colon_secretion_gene = {''}; % should become an input questions
        [model] = addReactionsHH(model, {['Colon_',upper(Metabolite),'t[bpC]']},strcat('Transport of ',space,MetaboliteName,' into the portal blood (Colon)'), {[substr,' <==> ',prod]}, {Colon_secretion_gene}, {'Transport, extracellular'},couplingFactor);
    end
    % excretion into bp
    % 'Colon_EX_glc_D[bpC]_[bp]'	'Exchange of D-Glucose (Colon)'	'Colon_glc_D[bpC]  -> glc_D[bp] '	'Transport, biofluid'
    substr = ['Colon_',Metabolite,'[bpC]'];
    prod = [Metabolite,'[bp]'];
    [model] = addReactionsHH(model, {['Colon_EX_',Metabolite,'[bpC]_[bp]']},strcat('Transport of ',space,MetaboliteName,' from Extracellular Space to portal vein (Colon)'), {[substr,' --> ',prod]}, {''}, {'Transport, biofluid'},couplingFactor);
end
%% transport into liver from bp and excretion into blood
if (isfield(param,'SIabsorp')  && ~isempty(param.SIabsorp) && contains(lower(param.SIabsorp),'y'))...
        || (isfield(param,'LIabsorp')  && ~isempty(param.LIabsorp) && contains(lower(param.LIabsorp),'y')) ...
        || (isfield(param,'PeripheralOrgan')  && any(ismember(param.PeripheralOrgan,'Gall'))) ...
        || (isfield(param,'PeripheralOrgan')  && any(ismember(param.PeripheralOrgan,'sIEC'))) ...
        || (isfield(param,'PeripheralOrgan')  && any(ismember(param.PeripheralOrgan,'Colon')))  ...
        || (isfield(param,'PeripheralOrgan')  && any(ismember(param.PeripheralOrgan,'Pancreas')))
    
    % 'Liver_EX_glc_D[bpL]_[bp]'	'Exchange of D-Glucose (Liver)'	'glc_D[bp]  -> Liver_glc_D[bpL] '		'Transport, biofluid'
    %
        substr = [Metabolite,'[bp]'];
        prod = ['Liver_',Metabolite,'[bpL]'];
        [model] = addReactionsHH(model, {['Liver_EX_',Metabolite,'[bpL]_[bp]']},strcat('Exchange of',space,MetaboliteName,' (Liver)'), {[substr,' --> ',prod]}, {''}, {'Transport, biofluid'},couplingFactor);

    % transport from the bp
    % e.g., note that gene is given here
    % 'Liver_GLCt1r[bpL]'	'Glucose Transport (Uniport) (Liver)'	'Liver_glc_D[bpL]  <=> Liver_glc_D[c] ''66035.1 or 6517.1 or 56606.1 or 154091.1 or 6513.1 or 11182.1 or 29988.1 or 6514.1 or 155184.1 or 56606.2 or 144195.1 or 6515.1 ...'	'Transport, extracellular'
    % check if reaction is specified in reaction structure
    if exist('reaction','var') && isfield(reaction,'Liver')
        organ = 'Liver';
        % find all relevant reactions
        e = reaction.(organ).rxns.abbr(contains(reaction.(organ).rxns.formula,'[bpL]'));
        c = reaction.(organ).rxns.abbr(contains(reaction.(organ).rxns.formula,'[c]'));
        rxnT = intersect(e,c);
        if ~isempty(rxnT)
            model = addOrganRxn(model,reaction, organ,rxnT, couplingFactor,'[bpL]');
        else
            % default reaction
            substr = ['Liver_',Metabolite,'[bpL]'];
            prod = ['Liver_',Metabolite,'[c]'];
            Liver_uptake_gene = {''}; % should become an input questions
            [model] = addReactionsHH(model, {['Liver_',upper(Metabolite),'t[bpL]']},strcat('Transport of ',space,MetaboliteName,' from Extracellular Space (portal vein) to Cytosol of Liver'), {[substr,' <==> ',prod]}, {Liver_uptake_gene}, {'Transport, extracellular'},couplingFactor);
            
        end
    else
        % default reaction
        substr = ['Liver_',Metabolite,'[bpL]'];
        prod = ['Liver_',Metabolite,'[c]'];
        Liver_uptake_gene = {''}; % should become an input questions
        [model] = addReactionsHH(model, {['Liver_',upper(Metabolite),'t[bpL]']},strcat('Transport of ',space,MetaboliteName,' from Extracellular Space (portal vein) to Cytosol of Liver'), {[substr,' <==> ',prod]}, {Liver_uptake_gene}, {'Transport, extracellular'},couplingFactor);
    end
    % check if reaction is specified in reaction structure
    if exist('reaction','var') && isfield(reaction,'Liver')
        organ = 'Liver';
        % find all relevant reactions
        e = reaction.(organ).rxns.abbr(contains(reaction.(organ).rxns.formula,'[e]'));
        c = reaction.(organ).rxns.abbr(contains(reaction.(organ).rxns.formula,'[c]'));
        rxnT = intersect(e,c);
        if ~isempty(rxnT)
            % if this is an abc transporter - also enable import
            % with default reaction
            if (contains(rxnT,'ABCt'))
                % check that indeed an abc efflux
                for r = 1: length(rxnT)
                    x = find(ismember(reaction.(organ).rxns.abbr,rxnT{r}));
                    if ~isempty(contains(   reaction.(organ).rxns.formula(x),'atp[c]'))
                        % add to abbr [bc]
                        %        reaction.(organ).rxns.abbr(x)=strcat( reaction.(organ).rxns.abbr(x),'[bc]');
                        reaction.(organ).rxns.abbr(x)=strcat( reaction.(organ).rxns.abbr(x));
                        % default reaction
                        substr = [organ,'_',Metabolite,'[c]'];
                        prod= [organ,'_',Metabolite,'[e]'];
                        Organ_secretion_gene = {''}; % should become an input questions
                        %                [model] = addReactionsHH(model, {[organ,'_',upper(Metabolite),'t[bc]']},strcat('Transport of ',space,MetaboliteName,' from Cytosol to Extracellular Space (blood) (',organ,')'), {[substr,' <==> ',prod]}, {Organ_secretion_gene}, {'Transport, extracellular'},couplingFactor);
                        [model] = addReactionsHH(model, {[organ,'_',upper(Metabolite),'t']},strcat('Transport of ',space,MetaboliteName,' from Cytosol to Extracellular Space (blood) (',organ,')'), {[substr,' <==> ',prod]}, {Organ_secretion_gene}, {'Transport, extracellular'},couplingFactor);
                    end
                end
            end
            %    model = addOrganRxn(model,reaction, organ,rxnT, couplingFactor,'[bc]');
            model = addOrganRxn(model,reaction, organ,rxnT, couplingFactor);
        else % default reaction
            % This is a default reaction to be added if no info is provided
            substr = ['Liver_',Metabolite,'[c]'];
            prod= ['Liver_',Metabolite,'[e]'];
            Liver_secretion_gene = {''}; % should become an input questions
            %      [model] = addReactionsHH(model, {['Liver_',upper(Metabolite),'t[bc]']},strcat('Transport of ',space,MetaboliteName,' from Cytosol to Extracellular Space (blood) of Liver'), {[substr,' --> ',prod]}, {Liver_uptake_gene}, {'Transport, extracellular'},couplingFactor);
            [model] = addReactionsHH(model, {['Liver_',upper(Metabolite),'t']},strcat('Transport of ',space,MetaboliteName,' from Cytosol to Extracellular Space (blood) of Liver'), {[substr,' --> ',prod]}, {Liver_uptake_gene}, {'Transport, extracellular'},couplingFactor);
        end
    else
        % This is a default reaction to be added if no info is provided
        substr = ['Liver_',Metabolite,'[c]'];
        prod= ['Liver_',Metabolite,'[e]'];
        Liver_secretion_gene = {''}; % should become an input questions
        %    [model] = addReactionsHH(model, {['Liver_',upper(Metabolite),'t[bc]']},strcat('Transport of ',space,MetaboliteName,' from Cytosol to Extracellular Space (blood) of Liver'), {[substr,' --> ',prod]}, {Liver_uptake_gene}, {'Transport, extracellular'},couplingFactor);
        [model] = addReactionsHH(model, {['Liver_',upper(Metabolite),'t']},strcat('Transport of ',space,MetaboliteName,' from Cytosol to Extracellular Space (blood) of Liver'), {[substr,' --> ',prod]}, {Liver_uptake_gene}, {'Transport, extracellular'},couplingFactor);
    end
    % 'Liver_EX_glc_D(e)_[bc]'	'Liver_EX_glc_D(e)_[bc]'	'Liver_glc_D[e]  -> glc_D[bc] ' 'Transport, biofluid'
    substr = ['Liver_',Metabolite,'[e]'];
    prod = [Metabolite,'[bc]'];
    [model] = addReactionsHH(model, {['Liver_EX_',Metabolite,'(e)_[bc]']},strcat('Exchange of',space,MetaboliteName,' (Liver)'), {[substr,' <==> ',prod]}, {''}, {'Transport, biofluid'},couplingFactor);
    
end

%% transport from Liver to bd
if (isfield(param,'BileDuct')  && ~isempty(param.BileDuct) && contains(lower(param.BileDuct),'y')) %|| (~isempty(param.LIabsorp) && contains(lower(param.LIabsorp),'y'))
    
    % 'Liver_C10CRNe[bdL]'	'Transport of Decanoyl Carnitine from Extra Cellular Space (Liver, bile duct)'	'Liver_c10crn[c]  -> Liver_c10crn[bdL] '	'Fatty acid oxidation'
    substr = ['Liver_',Metabolite,'[c]'];
    prod = ['Liver_',Metabolite,'[bdL]'];
    [model] = addReactionsHH(model, {['Liver_EX_',Metabolite,'[bdL]']},strcat('Exchange of',space,MetaboliteName,' (Liver, bile duct)'), {[substr,' --> ',prod]}, {''}, {'Transport, biofluid'},couplingFactor);
    
    % 'Liver_EX_c10crn[bdL]_[bd]'	'Exchange of Decanoyl Carnitine (Liver)'	'Liver_c10crn[bdL]  -> c10crn[bd] '		'Transport, biofluid'
    substr = ['Liver_',Metabolite,'[bdL]'];
    prod = [Metabolite,'[bd]'];
    [model] = addReactionsHH(model, {['Liver_EX_',Metabolite,'[bdL]_[bd]']},strcat('Exchange of',space,MetaboliteName,' (Liver, bile duct)'), {[substr,' --> ',prod]}, {''}, {'Transport, biofluid'},couplingFactor);
    
    %% transport from bd to LuSI
    % 'BileDuct_EX_c10crn[bd]_[luSI]'	'Exchange of Decanoyl Carnitine (BileDuct)'	'c10crn[bd]  -> c10crn[luSI] '	'Transport, biofluid'
    
    substr = [Metabolite,'[bd]'];
    prod = [Metabolite,'[luSI]'];
    [model] = addReactionsHH(model, {['BileDuct_EX_',Metabolite,'[bd]_[luSI]']},strcat('Exchange of',space,MetaboliteName,' (Bile duct)'), {[substr,' --> ',prod]}, {''}, {'Transport, biofluid'},couplingFactor);
    
    % 'SI_EX_glc_D[luSI]_[luLI]'	'Exchange of D-Glucose (SI)'	'glc_D[luSI]  -> glc_D[luLI] '	'Transport, biofluid'
    substr = [Metabolite,'[luSI]'];
    prod = [Metabolite,'[luLI]'];
    [model] = addReactionsHH(model, {['SI_EX_',Metabolite,'[luSI]_[luLI]']},strcat('Exchange of',space,MetaboliteName,' (SI)'), {[substr,' --> ',prod]}, {''}, {'Transport, biofluid'},couplingFactor);
    
    % 'LI_EX_glc_D[luLI]_[fe]'	'Exchange of D-Glucose (LI)'	'glc_D[luLI]  -> glc_D[fe] '	'Transport, biofluid'
    substr = [Metabolite,'[luLI]'];
    prod = [Metabolite,'[fe]'];
    [model] = addReactionsHH(model, {['LI_EX_',Metabolite,'[luLI]_[fe]']},strcat('Exchange of',space,MetaboliteName,' (LI)'), {[substr,' --> ',prod]}, {''}, {'Transport, biofluid'},couplingFactor);
    
    % 'Excretion_EX_glc_D[fe]'	'Exchange of D-Glucose (Excretion)'	'glc_D[fe]  -> '	'Exchange/demand reaction'
    substr = [Metabolite,'[fe]'];
    [model] = addReactionsHH(model, {['Excretion_EX_',Metabolite,'[fe]']},strcat('Exchange of',space,MetaboliteName,' (Excretion)'), {[substr,' --> ']}, {''}, {'Exchange/demand reaction'},couplingFactor);
    
end

%% metabolite goes straight into blood stream - IV
if (isfield(param,'Blood') && ~isempty(param.Blood) && contains(lower(param.Blood),'y'))
    % 1. 'Diet_EX_glc_D[d]'	'Exchange of D-Glucose (Diet)'	'glc_D[d] -->' 'Exchange/demand reaction'
    substr = [Metabolite,'[ext]']; % new compartment [ext] external
    [model] = addReactionsHH(model, {['EX_',Metabolite,'[ext]']},strcat('Exchange of',space,MetaboliteName,' (External)'), {[substr,' -->']}, {''}, {'Exchange/demand reaction'},couplingFactor);
    % 'GI_EX_glc_D[lu]_[d]'	'Exchange of D-Glucose (GI)'	'glc_D[d]  -> glc_D[lu] '   'Transport, biofluid'
    substr = [Metabolite,'[ext]'];
    prod = [Metabolite,'[bc]'];
    [model] = addReactionsHH(model, {['BC_EX_',Metabolite,'[ext]_[bc]']},strcat('Exchange of',space,MetaboliteName,' (External)'), {[substr,' --> ',prod]}, {''}, {'Transport, biofluid'},couplingFactor);
end
%% metabolite goes straight into blood stream - transdermal
if (isfield(param,'Transdermal') && ~isempty(param.Transdermal) && contains(lower(param.Transdermal),'y'))
    % 1. 'Diet_EX_glc_D[d]'	'Exchange of D-Glucose (Diet)'	'glc_D[d] -->' 'Exchange/demand reaction'
    substr = [Metabolite,'[ext]']; % new compartment [ext] external
    [model] = addReactionsHH(model, {['EX_',Metabolite,'[ext]']},strcat('Exchange of',space,MetaboliteName,' (External)'), {[substr,' -->']}, {''}, {'Exchange/demand reaction'},couplingFactor);
    % 'GI_EX_glc_D[lu]_[d]'	'Exchange of D-Glucose (GI)'	'glc_D[d]  -> glc_D[lu] '   'Transport, biofluid'
    substr = [Metabolite,'[ext]'];
    prod = [Metabolite,'[bc]'];
    [model] = addReactionsHH(model, {['BC_EX_',Metabolite,'[ext]_[bc]']},strcat('Exchange of',space,MetaboliteName,' (External)'), {[substr,' --> ',prod]}, {''}, {'Transport, biofluid'},couplingFactor);
end

%% % topical: inhalation (TI)
if (isfield(param,'Inhalation') && ~isempty(param.Inhalation) && contains(lower(param.Inhalation),'y'))
    %EX_h2o[a]	h2o[a] 	->
    substr = [Metabolite,'[a]']; % new compartment [ext] external
    [model] = addReactionsHH(model, {['EX_',Metabolite,'[a]']},strcat('Exchange of',space,MetaboliteName,' (Air)'), {[substr,' <==>']}, {''}, {'Exchange/demand reaction'},couplingFactor);
    %{'Lung_EX_h2o(a)_[aL]'}     {'Lung_h2o[aL]  <=> h2o[a] '}
    substr = ['Lung_',Metabolite,'[aL]'];
    prod = [Metabolite,'[a]'];
    [model] = addReactionsHH(model, {['Lung_EX_',Metabolite,'(a)_[aL]']},strcat('Exchange of',space,MetaboliteName,' (Air)'), {[substr,' <==> ',prod]}, {''}, {'Transport, biofluid'},couplingFactor);
    %  {'Lung_H2Ot[aL]'      } {'Lung_h2o[aL]  <=> Lung_h2o[c] '}
    substr = ['Lung_',Metabolite,'[aL]'];
    prod = ['Lung_',Metabolite,'[c]'];
    [model] = addReactionsHH(model, {['Lung_',upper(Metabolite),'t[aL]']},strcat('Transport of',space,MetaboliteName,' (from aL to c)'), {[substr,' <==> ',prod]}, {''}, {'Transport, biofluid'},couplingFactor);
    %   {'Lung_H2Ot'                 }'Lung_h2o[e]  <=> Lung_h2o[c] '
    substr = ['Lung_',Metabolite,'[e]'];
    prod = ['Lung_',Metabolite,'[c]'];
    %   [model] = addReactionsHH(model, {['Lung_',upper(Metabolite),'t[bc]']},strcat('Transport of',space,MetaboliteName,' (from c to e)'), {[substr,' <==> ',prod]}, {''}, {'Transport, biofluid'},couplingFactor);
    [model] = addReactionsHH(model, {['Lung_',upper(Metabolite),'t']},strcat('Transport of',space,MetaboliteName,' (from c to e)'), {[substr,' <==> ',prod]}, {''}, {'Transport, biofluid'},couplingFactor);
    %      {'Lung_EX_h2o(e)_[bc]'}}{'Lung_h2o[e]  <=> h2o[bc] '   substr = [Metabolite,'[aL]'];
    substr = ['Lung_',Metabolite,'[e]'];
    prod = [Metabolite,'[bc]'];
    [model] = addReactionsHH(model, {['Lung_EX_',Metabolite,'(e)_[bc]']},strcat('Transport of',space,MetaboliteName,' (from e to bc)'), {[substr,' <==> ',prod]}, {''}, {'Transport, biofluid'},couplingFactor);
end

%% metabolite goes straight into blood stream from Diet 
if (isfield(param,'Sublingual') && ~isempty(param.Sublingual) && contains(lower(param.Sublingual),'y'))
    % 1. 'Diet_EX_glc_D[d]'	'Exchange of D-Glucose (Diet)'	'glc_D[d] -->' 'Exchange/demand reaction'
    substr = [Metabolite,'[d]']; % new compartment [ext] external
    [model] = addReactionsHH(model, {['Diet_EX_',Metabolite,'[d]']},strcat('Exchange of',space,MetaboliteName,' (d)'), {[substr,' -->']}, {''}, {'Exchange/demand reaction'},couplingFactor);
    % 'GI_EX_glc_D[lu]_[d]'	'Exchange of D-Glucose (GI)'	'glc_D[d]  -> glc_D[lu] '   'Transport, biofluid'
    substr = [Metabolite,'[d]'];
    prod = [Metabolite,'[bc]'];
    [model] = addReactionsHH(model, {['BC_EX_',Metabolite,'[d]_[bc]']},strcat('Exchange of',space,MetaboliteName,' (Diet to blood)'), {[substr,' --> ',prod]}, {''}, {'Transport, biofluid'},couplingFactor);
end
%% breathing out of volatives (produced by microbes)
if (isfield(param,'Breath') && ~isempty(param.Breath) && contains(lower(param.Breath),'y'))
    % from LI to SI
    % 'SI_EX_glc_D[luSI]_[luLI]'	'Exchange of D-Glucose (SI)'	'glc_D[luSI]  -> glc_D[luLI] '	'Transport, biofluid'
    substr = [Metabolite,'[luSI]'];
    prod = [Metabolite,'[luLI]'];
    % check if reaction already exist -- if so, set to reversible
    clear rxnA
    rxnA = find(contains(model.rxns,['SI_EX_',Metabolite,'[luSI]_[luLI]']));
    if isempty(rxnA)
        [model] = addReactionsHH(model, {['SI_EX_',Metabolite,'[luSI]_[luLI]']},strcat('Exchange of',space,MetaboliteName,' (LI to SI)'), {[substr,' <==> ',prod]}, {''}, {'Transport, biofluid'},couplingFactor);
    else
        model.lb(rxnA) = -1*model.ub(rxnA);
    end
    % from SI to lu
    % 'SI_EX_glc_D[lu]_[luSI]'	'Exchange of D-Glucose (SI)'	'glc_D[lu]  -> glc_D[luSI] ' 'Transport, biofluid'
    clear rxnA
    substr = [Metabolite,'[lu]'];
    prod = [Metabolite,'[luSI]'];
    rxnA = find(contains(model.rxns,['SI_EX_',Metabolite,'[lu]_[luSI]']));
    if isempty(rxnA)
        [model] = addReactionsHH(model, {['SI_EX_',Metabolite,'[lu]_[luSI]']},strcat('Exchange of',space,MetaboliteName,' (SI to LU)'), {[substr,' <==> ',prod]}, {''}, {'Transport, biofluid'},couplingFactor);
    else
        model.lb(rxnA) = -1*model.ub(rxnA);
    end
    %    from lu to a
    substr = [Metabolite,'[a]'];
    prod = [Metabolite,'[lu]'];
    [model] = addReactionsHH(model, {['GI_EX_',Metabolite,'[lu]_[a]']},strcat('Exchange of',space,MetaboliteName,' (GI to)'), {[substr,' <==> ',prod]}, {''}, {'Transport, biofluid'},couplingFactor);
    
    %EX_h2o[a]	h2o[a] 	->
    substr = [Metabolite,'[a]']; % new compartment [ext] external
    [model] = addReactionsHH(model, {['EX_',Metabolite,'[a]']},strcat('Exchange of',space,MetaboliteName,' (Air)'), {[substr,' -->']}, {''}, {'Exchange/demand reaction'},couplingFactor);
end



%% kidney uptake
if (isfield(param,'MicrobiomeOnly') && ~isempty(param.MicrobiomeOnly) && contains(lower(param.MicrobiomeOnly),'y')) ...
        || (~isfield(param,'PeripheralOrgan') || isempty(param.PeripheralOrgan))
    % no kidney reactions
else
    % 'Kidney_EX_glc_D(e)_[bc]'	'Exchange of D-Glucose (frombloodto[e]'	'glc_D[e]  <=> Kidney_glc_D[bc] '	'Transport, biofluid'
    substr = [Metabolite,'[bc]'];
    prod = ['Kidney_',Metabolite,'[e]'];
    
    [model] = addReactionsHH(model, {['Kidney_EX_',Metabolite,'(e)_[bc]']},strcat('Exchange of',space,MetaboliteName,' (from blood to [e])'), {[substr,' <==> ',prod]}, {''}, {'Transport, biofluid'},couplingFactor);
    
    %% reabsorption
    % 'Kidney_EX_glc_D[bcK]_[bc]'	'Exchange of D-Glucose (Kidney)'	'Kidney_glc_D[bcK]  -> glc_D[bc] '	0	1000000	{''}	'Transport, biofluid'
    substr = ['Kidney_',Metabolite,'[bcK]'];
    prod = [Metabolite,'[bc]'];
    [model] = addReactionsHH(model, {['Kidney_EX_',Metabolite,'[bcK]_[bc]']},strcat('Exchange of',space,MetaboliteName,' (Kidney)'), {[substr,' --> ',prod]}, {''}, {'Transport, biofluid'},couplingFactor);
    
    % 'Kidney_EX_11docrtsl(e)_[bcK]'	'Kidney_11docrtsl(e)_[bcK]'	'Kidney_11docrtsl[e]  -> Kidney_11docrtsl[bcK] '	0	1000000	{''}	'Transport, passive reabsorption'
    substr = ['Kidney_',Metabolite,'[e]'];
    prod = ['Kidney_',Metabolite,'[bcK]'];
    [model] = addReactionsHH(model, {['Kidney_EX_',Metabolite,'(e)_[bcK]']},strcat('Exchange of',space,MetaboliteName,' (Kidney)'), {[substr,' --> ',prod]}, {''}, {'Transport, biofluid'},couplingFactor);
    
    %% urine excretion
    % 'Kidney_EX_glc_D(e)_[u]'	'Exchange of D-Glucose (from [e] to blood)'	'Kidney_glc_D[e]  -> glc_D[u] '		'Transport, biofluid'
    substr = ['Kidney_',Metabolite,'[e]'];
    prod = [Metabolite,'[u]'];
    [model] = addReactionsHH(model, {['Kidney_EX_',Metabolite,'(e)_[u]']},strcat('Exchange of',space,MetaboliteName,' (from [e] to urine)'), {[substr,' --> ',prod]}, {''}, {'Transport, biofluid'},couplingFactor);
    
    % 'EX_glc_D[u]'	'EX_glc_D[u]'	'glc_D[u]  -> ' 	{''}	'Exchange/demand reaction'
    substr = [Metabolite,'[u]'];
    prod = [''];
    [model] = addReactionsHH(model, {['EX_',Metabolite,'[u]']},strcat('Excretion of',space,MetaboliteName,' (urine)'), {[substr,' --> ',prod]}, {''}, {'Exchange/demand reaction'},couplingFactor);
end
% peripheral organ
if (isfield(param,'PeripheralOrgan') && ~isempty(param.PeripheralOrgan))
    for i = 1 :length(param.PeripheralOrgan)
        if strcmp(param.PeripheralOrgan{i},'Pancreas')
            if exist('reaction','var') && isfield(reaction,param.PeripheralOrgan{i})
                organ = param.PeripheralOrgan{i};
                % transport into the Pancreas
                % find all relevant reactions
                e = reaction.(organ).rxns.abbr(contains(reaction.(organ).rxns.formula,'[e]'));
                c = reaction.(organ).rxns.abbr(contains(reaction.(organ).rxns.formula,'[c]'));
                rxnT = intersect(e,c);
                if ~isempty(rxnT)
                    % if this is an abc transporter - also enable import
                    % with default reaction
                    if (contains(rxnT,'ABCt'))
                        % check that indeed an abc efflux
                        for r = 1: length(rxnT)
                            x = find(ismember(reaction.(organ).rxns.abbr,rxnT{r}));
                            if ~isempty(contains(   reaction.(organ).rxns.formula(x),'atp[c]'))
                                % add to abbr [bc]
                                %   reaction.(organ).rxns.abbr(x)=strcat( reaction.(organ).rxns.abbr(x),'[bc]');
                                reaction.(organ).rxns.abbr(x)=strcat( reaction.(organ).rxns.abbr(x));
                                % default reaction
                                substr = [param.PeripheralOrgan{i},'_',Metabolite,'[c]'];
                                prod= [param.PeripheralOrgan{i},'_',Metabolite,'[e]'];
                                Organ_secretion_gene = {''}; % should become an input questions
                                %      [model] = addReactionsHH(model, {[param.PeripheralOrgan{i},'_',upper(Metabolite),'t[bc]']},strcat('Transport of ',space,MetaboliteName,' from Cytosol to Extracellular Space (blood) (',param.PeripheralOrgan{i},')'), {[substr,' <==> ',prod]}, {Organ_secretion_gene}, {'Transport, extracellular'},couplingFactor);
                                [model] = addReactionsHH(model, {[param.PeripheralOrgan{i},'_',upper(Metabolite),'t']},strcat('Transport of ',space,MetaboliteName,' from Cytosol to Extracellular Space (blood) (',param.PeripheralOrgan{i},')'), {[substr,' <==> ',prod]}, {Organ_secretion_gene}, {'Transport, extracellular'},couplingFactor);
                            end
                        end
                    end
                    %        model = addOrganRxn(model,reaction, organ,rxnT, couplingFactor,'[bc]');
                    model = addOrganRxn(model,reaction, organ,rxnT, couplingFactor);
                else
                    % default reaction
                    substr = [param.PeripheralOrgan{i},'_',Metabolite,'[c]'];
                    prod= [param.PeripheralOrgan{i},'_',Metabolite,'[e]'];
                    
                    Organ_secretion_gene = {''}; % should become an input questions
                    %          [model] = addReactionsHH(model, {[param.PeripheralOrgan{i},'_',upper(Metabolite),'t[bc]']},strcat('Transport of ',space,MetaboliteName,' from Cytosol to Extracellular Space (blood) (',param.PeripheralOrgan{i},')'), {[substr,' <==> ',prod]}, {Organ_secretion_gene}, {'Transport, extracellular'},couplingFactor);
                    [model] = addReactionsHH(model, {[param.PeripheralOrgan{i},'_',upper(Metabolite),'t']},strcat('Transport of ',space,MetaboliteName,' from Cytosol to Extracellular Space (blood) (',param.PeripheralOrgan{i},')'), {[substr,' <==> ',prod]}, {Organ_secretion_gene}, {'Transport, extracellular'},couplingFactor);
                end
            else % default reaction
                substr = [param.PeripheralOrgan{i},'_',Metabolite,'[c]'];
                prod= [param.PeripheralOrgan{i},'_',Metabolite,'[e]'];
                
                Organ_secretion_gene = {''}; % should become an input questions
                %   [model] = addReactionsHH(model, {[param.PeripheralOrgan{i},'_',upper(Metabolite),'t[bc]']},strcat('Transport of ',space,MetaboliteName,' from Cytosol to Extracellular Space (blood) (',param.PeripheralOrgan{i},')'), {[substr,' <==> ',prod]}, {Organ_secretion_gene}, {'Transport, extracellular'},couplingFactor);
                [model] = addReactionsHH(model, {[param.PeripheralOrgan{i},'_',upper(Metabolite),'t']},strcat('Transport of ',space,MetaboliteName,' from Cytosol to Extracellular Space (blood) (',param.PeripheralOrgan{i},')'), {[substr,' <==> ',prod]}, {Organ_secretion_gene}, {'Transport, extracellular'},couplingFactor);
            end
            % 'Liver_EX_glc_D(e)_[bc]'	'Liver_EX_glc_D(e)_[bc]'	'Liver_glc_D[e]  -> glc_D[bc] ' 'Transport, biofluid'
            substr = [param.PeripheralOrgan{i},'_',Metabolite,'[e]'];
            prod = [Metabolite,'[bc]'];
            [model] = addReactionsHH(model, {[param.PeripheralOrgan{i},'_EX_',Metabolite,'(e)_[bc]']},strcat('Exchange of',space,MetaboliteName,' (',param.PeripheralOrgan{i},')'), {[substr,' <==> ',prod]}, {''}, {'Transport, biofluid'},couplingFactor);
            model.ub(ismember(model.rxns,[param.PeripheralOrgan{i},'_EX_',Metabolite,'(e)_[bc]'])) = 0; % uptake only
            % secretion is into [bp]
            if exist('reaction','var') && isfield(reaction,'Pancreas')
                % find all relevant reactions
                e = reaction.(organ).rxns.abbr(contains(reaction.(organ).rxns.formula,'[bpP]'));
                c = reaction.(organ).rxns.abbr(contains(reaction.(organ).rxns.formula,'[c]'));
                rxnT = intersect(e,c);
                if ~isempty(rxnT)
                    model = addOrganRxn(model,reaction, organ,rxnT, couplingFactor,'[bpP]');
                else
                    % default reaction
                    substr = ['Pancreas_',Metabolite,'[c]'];
                    prod = ['Pancreas_',Metabolite,'[bpP]'];
                    SIEC_secretion_gene = {''}; % should become an input questions
                    [model] = addReactionsHH(model, {['Pancreas_',upper(Metabolite),'t[bpP]']},strcat('Transport of ',space,MetaboliteName,' into the portal blood (bpP)'), {[substr,' --> ',prod]}, {SIEC_secretion_gene}, {'Transport, extracellular'},couplingFactor);
                end
            else
                % default reaction
                substr = ['Pancreas_',Metabolite,'[c]'];
                prod = ['Pancreas_',Metabolite,'[bpP]'];
                SIEC_secretion_gene = {''}; % should become an input questions
                [model] = addReactionsHH(model, {['Pancreas_',upper(Metabolite),'t[bpP]']},strcat('Transport of ',space,MetaboliteName,' into the portal blood (bpP)'), {[substr,' --> ',prod]}, {SIEC_secretion_gene}, {'Transport, extracellular'},couplingFactor);
            end
            % excretion into bp
            % 'sIEC_EX_glc_D[bpI]_[bp]'	'Exchange of D-Glucose (sIEC)'	'sIEC_glc_D[bpI]  -> glc_D[bp] '	'Transport, biofluid'
            substr = ['Pancreas_',Metabolite,'[bpP]'];
            prod = [Metabolite,'[bp]'];
            [model] = addReactionsHH(model, {['Pancreas_EX_',Metabolite,'[bpP]_[bp]']},strcat('Transport of ',space,MetaboliteName,' from Extracellular Space to portal vein (Pancreas)'), {[substr,' --> ',prod]}, {''}, {'Transport, biofluid'},couplingFactor);
            
        elseif strcmp(param.PeripheralOrgan{i},'Brain') || strcmp(param.PeripheralOrgan{i},'Scord')
            % Brain_TCMNt[csf]: Brain_tcmn[c] 	->	Brain_tcmn[e]
            if exist('reaction','var') && isfield(reaction,param.PeripheralOrgan{i})
                organ = param.PeripheralOrgan{i};
                % find all relevant reactions
                e = reaction.(organ).rxns.abbr(contains(reaction.(organ).rxns.formula,'[e]'));
                c = reaction.(organ).rxns.abbr(contains(reaction.(organ).rxns.formula,'[c]'));
                rxnT = intersect(e,c);
                if ~isempty(rxnT)
                    % if this is an abc transporter - also enable import
                    % with default reaction
                    if (contains(rxnT,'ABCt'))
                        % check that indeed an abc efflux
                        for r = 1: length(rxnT)
                            x = find(ismember(reaction.(organ).rxns.abbr,rxnT{r}));
                            if ~isempty(contains(   reaction.(organ).rxns.formula(x),'atp[c]'))
                                % add to abbr [bc]
                                %   reaction.(organ).rxns.abbr(x)=strcat( reaction.(organ).rxns.abbr(x),'[csf]');
                                reaction.(organ).rxns.abbr(x)=strcat( reaction.(organ).rxns.abbr(x));
                                % default reaction
                                substr = [param.PeripheralOrgan{i},'_',Metabolite,'[c]'];
                                prod= [param.PeripheralOrgan{i},'_',Metabolite,'[e]'];
                                Organ_secretion_gene = {''}; % should become an input questions
                                %          [model] = addReactionsHH(model, {[param.PeripheralOrgan{i},'_',upper(Metabolite),'t[csf]']},strcat('Transport of ',space,MetaboliteName,' from Cytosol to Extracellular Space (csf) (',param.PeripheralOrgan{i},')'), {[substr,' <==> ',prod]}, {Organ_secretion_gene}, {'Transport, extracellular'},couplingFactor);
                                [model] = addReactionsHH(model, {[param.PeripheralOrgan{i},'_',upper(Metabolite),'t']},strcat('Transport of ',space,MetaboliteName,' from Cytosol to Extracellular Space (csf) (',param.PeripheralOrgan{i},')'), {[substr,' <==> ',prod]}, {Organ_secretion_gene}, {'Transport, extracellular'},couplingFactor);
                            end
                        end
                    end
                    
                    %     model = addOrganRxn(model,reaction, organ,rxnT, couplingFactor,'[csf]');
                    model = addOrganRxn(model,reaction, organ,rxnT, couplingFactor);
                else
                    % default reaction
                    substr = [param.PeripheralOrgan{i},'_',Metabolite,'[c]'];
                    prod= [param.PeripheralOrgan{i},'_',Metabolite,'[e]'];
                    Organ_secretion_gene = {''}; % should become an input questions
                    %       [model] = addReactionsHH(model, {[param.PeripheralOrgan{i},'_',upper(Metabolite),'t[csf]']},strcat('Transport of ',space,MetaboliteName,' from Cytosol to Extracellular Space (csf) (',param.PeripheralOrgan{i},')'), {[substr,' <==> ',prod]}, {Organ_secretion_gene}, {'Transport, extracellular'},couplingFactor);
                    [model] = addReactionsHH(model, {[param.PeripheralOrgan{i},'_',upper(Metabolite),'t']},strcat('Transport of ',space,MetaboliteName,' from Cytosol to Extracellular Space (csf) (',param.PeripheralOrgan{i},')'), {[substr,' <==> ',prod]}, {Organ_secretion_gene}, {'Transport, extracellular'},couplingFactor);
                    
                end
            else
                % default reaction
                substr = [param.PeripheralOrgan{i},'_',Metabolite,'[c]'];
                prod= [param.PeripheralOrgan{i},'_',Metabolite,'[e]'];
                Organ_secretion_gene = {''}; % should become an input questions
                %    if length(find(ismember(model.rxns, [param.PeripheralOrgan{i},'_',upper(Metabolite),'t[csf]'])))>0 % reaction does not exist in model
                %         [model] = addReactionsHH(model, {[param.PeripheralOrgan{i},'_',upper(Metabolite),'t[csf]']},strcat('Transport of ',space,MetaboliteName,' from Cytosol to Extracellular Space (csf) (',param.PeripheralOrgan{i},')'), {[substr,' <==> ',prod]}, {Organ_secretion_gene}, {'Transport, extracellular'},couplingFactor);
               % if length(find(ismember(model.rxns, [param.PeripheralOrgan{i},'_',upper(Metabolite),'t'])))>0 % reaction does not exist in model
                    [model] = addReactionsHH(model, {[param.PeripheralOrgan{i},'_',upper(Metabolite),'t[csf]']},strcat('Transport of ',space,MetaboliteName,' from Cytosol to Extracellular Space (csf) (',param.PeripheralOrgan{i},')'), {[substr,' <==> ',prod]}, {Organ_secretion_gene}, {'Transport, extracellular'},couplingFactor);
               % end
            end
            % 'Brain_EX_glc_D(e)_[csf]'
            % 'Scord_EX_glc_D(e)_[csf]'
            substr = [param.PeripheralOrgan{i},'_',Metabolite,'[e]'];
            prod = [Metabolite,'[csf]'];
            [model] = addReactionsHH(model, {[param.PeripheralOrgan{i},'_EX_',Metabolite,'(e)_[csf]']},strcat('Exchange of',space,MetaboliteName,' (',param.PeripheralOrgan{i},')'), {[substr,' <==> ',prod]}, {''}, {'Transport, biofluid'},couplingFactor);
            
            % 'BBB_TYM[CSF]exp'
            substr = [Metabolite,'[csf]'];
            prod = [Metabolite,'[bc]'];
            [model] = addReactionsHH(model, {['BBB_',upper(Metabolite),'[CSF]exp']},strcat('Export of',space,MetaboliteName,' from csf to bc'), {[substr,' --> ',prod]}, {''}, {'Transport, biofluid'},couplingFactor);
            
            % 'BBB_TYM[CSF]upt'
            substr = [Metabolite,'[csf]'];
            prod = [Metabolite,'[bc]'];
            [model] = addReactionsHH(model, {['BBB_',upper(Metabolite),'[CSF]upt']},strcat('Transport of',space,MetaboliteName,' across BBB'), {[substr,' <==> ',prod]}, {''}, {'Transport, biofluid'},couplingFactor);
            % reverse reaction only
            model.ub(contains(model.rxns,['BBB_',upper(Metabolite),'[CSF]upt']))=0;
        elseif strcmp(param.PeripheralOrgan{i},'Gall')
            % 'Gall_EX_glc_D[bdG]_[bd]'
            % gall can only uptake from [bc]
            if exist('reaction','var') && isfield(reaction,param.PeripheralOrgan{i})
                organ = param.PeripheralOrgan{i};
                % transport into the sIEC
                % e.g. note that gene is given here
                % 'sIEC_GLCSGLT1le[luI]'	'Transport of D-Glucose from Extracellular Space to Cytosol of Mucosal Cells in Small Intestine (sIEC)'	'sIEC_h2o[c] + 2 sIEC_na1[luI] + sIEC_atp[c] + sIEC_glc_D[luI]  -> sIEC_h[c] + 2 sIEC_na1[c] + sIEC_adp[c] + sIEC_pi[c] + sIEC_glc_D[c] '	'6523.1'	'Transport, extracellular'
                % !!!!!! also reaction may need to be specialised
                % find all relevant reactions
                e = reaction.(organ).rxns.abbr(contains(reaction.(organ).rxns.formula,'[e]'));
                c = reaction.(organ).rxns.abbr(contains(reaction.(organ).rxns.formula,'[c]'));
                rxnT = intersect(e,c);
                if ~isempty(rxnT)
                    % if this is an abc transporter - also enable import
                    % with default reaction
                    if (contains(rxnT,'ABCt'))
                        % check that indeed an abc efflux
                        for r = 1: length(rxnT)
                            x = find(ismember(reaction.(organ).rxns.abbr,rxnT{r}));
                            if ~isempty(contains(   reaction.(organ).rxns.formula(x),'atp[c]'))
                                % add to abbr [bc]
                                %     reaction.(organ).rxns.abbr(x)=strcat( reaction.(organ).rxns.abbr(x),'[bc]');
                                reaction.(organ).rxns.abbr(x)=strcat( reaction.(organ).rxns.abbr(x));
                                % default reaction
                                substr = [param.PeripheralOrgan{i},'_',Metabolite,'[c]'];
                                prod= [param.PeripheralOrgan{i},'_',Metabolite,'[e]'];
                                Organ_secretion_gene = {''}; % should become an input questions
                                %          [model] = addReactionsHH(model, {[param.PeripheralOrgan{i},'_',upper(Metabolite),'t[bc]']},strcat('Transport of ',space,MetaboliteName,' from Cytosol to Extracellular Space (blood) (',param.PeripheralOrgan{i},')'), {[substr,' <==> ',prod]}, {Organ_secretion_gene}, {'Transport, extracellular'},couplingFactor);
                                [model] = addReactionsHH(model, {[param.PeripheralOrgan{i},'_',upper(Metabolite),'t']},strcat('Transport of ',space,MetaboliteName,' from Cytosol to Extracellular Space (blood) (',param.PeripheralOrgan{i},')'), {[substr,' <==> ',prod]}, {Organ_secretion_gene}, {'Transport, extracellular'},couplingFactor);
                            end
                        end
                    end
                    %          model = addOrganRxn(model,reaction, organ,rxnT, couplingFactor,'[bc]');
                    model = addOrganRxn(model,reaction, organ,rxnT, couplingFactor);
                else
                    % default reaction
                    substr = [param.PeripheralOrgan{i},'_',Metabolite,'[c]'];
                    prod= [param.PeripheralOrgan{i},'_',Metabolite,'[e]'];
                    
                    Organ_secretion_gene = {''}; % should become an input questions
                    %           [model] = addReactionsHH(model, {[param.PeripheralOrgan{i},'_',upper(Metabolite),'t[bc]']},strcat('Transport of ',space,MetaboliteName,' from Cytosol to Extracellular Space (blood) (',param.PeripheralOrgan{i},')'), {[substr,' <==> ',prod]}, {Organ_secretion_gene}, {'Transport, extracellular'},couplingFactor);
                    [model] = addReactionsHH(model, {[param.PeripheralOrgan{i},'_',upper(Metabolite),'t']},strcat('Transport of ',space,MetaboliteName,' from Cytosol to Extracellular Space (blood) (',param.PeripheralOrgan{i},')'), {[substr,' <==> ',prod]}, {Organ_secretion_gene}, {'Transport, extracellular'},couplingFactor);
                end
            else % default reaction
                substr = [param.PeripheralOrgan{i},'_',Metabolite,'[c]'];
                prod= [param.PeripheralOrgan{i},'_',Metabolite,'[e]'];
                
                Organ_secretion_gene = {''}; % should become an input questions
                %        [model] = addReactionsHH(model, {[param.PeripheralOrgan{i},'_',upper(Metabolite),'t[bc]']},strcat('Transport of ',space,MetaboliteName,' from Cytosol to Extracellular Space (blood) (',param.PeripheralOrgan{i},')'), {[substr,' <==> ',prod]}, {Organ_secretion_gene}, {'Transport, extracellular'},couplingFactor);
                [model] = addReactionsHH(model, {[param.PeripheralOrgan{i},'_',upper(Metabolite),'t']},strcat('Transport of ',space,MetaboliteName,' from Cytosol to Extracellular Space (blood) (',param.PeripheralOrgan{i},')'), {[substr,' <==> ',prod]}, {Organ_secretion_gene}, {'Transport, extracellular'},couplingFactor);
            end
            % 'Liver_EX_glc_D(e)_[bc]'	'Liver_EX_glc_D(e)_[bc]'	'Liver_glc_D[e]  -> glc_D[bc] ' 'Transport, biofluid'
            substr = [param.PeripheralOrgan{i},'_',Metabolite,'[e]'];
            prod = [Metabolite,'[bc]'];
            [model] = addReactionsHH(model, {[param.PeripheralOrgan{i},'_EX_',Metabolite,'(e)_[bc]']},strcat('Exchange of',space,MetaboliteName,' (',param.PeripheralOrgan{i},')'), {[substr,' <==> ',prod]}, {''}, {'Transport, biofluid'},couplingFactor);
            model.ub(ismember(model.rxns,[param.PeripheralOrgan{i},'_EX_',Metabolite,'(e)_[bc]']))=0;
            % secretion is into [bp]
            if exist('reaction','var') && isfield(reaction,'Gall')
                % find all relevant reactions
                e = reaction.(organ).rxns.abbr(contains(reaction.(organ).rxns.formula,'[bpG]'));
                c = reaction.(organ).rxns.abbr(contains(reaction.(organ).rxns.formula,'[c]'));
                rxnT = intersect(e,c);
                if ~isempty(rxnT)
                    model = addOrganRxn(model,reaction, organ,rxnT, couplingFactor,'[bpG]');
                else
                    % default reaction
                    substr = ['Gall_',Metabolite,'[c]'];
                    prod = ['Gall_',Metabolite,'[bpG]'];
                    SIEC_secretion_gene = {''}; % should become an input questions
                    [model] = addReactionsHH(model, {['Gall_',upper(Metabolite),'t[bpG]']},strcat('Transport of ',space,MetaboliteName,' into the portal blood (bpG)'), {[substr,' --> ',prod]}, {SIEC_secretion_gene}, {'Transport, extracellular'},couplingFactor);
                end
            else
                % default reaction
                substr = ['Gall_',Metabolite,'[c]'];
                prod = ['Gall_',Metabolite,'[bpG]'];
                SIEC_secretion_gene = {''}; % should become an input questions
                [model] = addReactionsHH(model, {['Gall_',upper(Metabolite),'t[bpG]']},strcat('Transport of ',space,MetaboliteName,' into the portal blood (bpG)'), {[substr,' --> ',prod]}, {SIEC_secretion_gene}, {'Transport, extracellular'},couplingFactor);
            end
            % excretion into bp
            % 'sIEC_EX_glc_D[bpI]_[bp]'	'Exchange of D-Glucose (sIEC)'	'sIEC_glc_D[bpI]  -> glc_D[bp] '	'Transport, biofluid'
            substr = ['Gall_',Metabolite,'[bpG]'];
            prod = [Metabolite,'[bp]'];
            [model] = addReactionsHH(model, {['Gall_EX_',Metabolite,'[bpG]_[bp]']},strcat('Transport of ',space,MetaboliteName,' from Extracellular Space to portal vein (Gall)'), {[substr,' --> ',prod]}, {''}, {'Transport, biofluid'},couplingFactor);
            
        elseif strcmp(param.PeripheralOrgan{i},'Skin')
            if exist('reaction','var') && isfield(reaction,param.PeripheralOrgan{i})
                organ = param.PeripheralOrgan{i};
                % find all relevant reactions
                e = reaction.(organ).rxns.abbr(contains(reaction.(organ).rxns.formula,'[e]'));
                c = reaction.(organ).rxns.abbr(contains(reaction.(organ).rxns.formula,'[c]'));
                rxnT = intersect(e,c);
                if ~isempty(rxnT)
                    % if this is an abc transporter - also enable import
                    % with default reaction
                    if (contains(rxnT,'ABCt'))
                        % check that indeed an abc efflux
                        for r = 1: length(rxnT)
                            x = find(ismember(reaction.(organ).rxns.abbr,rxnT{r}));
                            if ~isempty(contains(   reaction.(organ).rxns.formula(x),'atp[c]'))
                                % add to abbr [bc]
                                %   reaction.(organ).rxns.abbr(x)=strcat( reaction.(organ).rxns.abbr(x),'[bc]');
                                reaction.(organ).rxns.abbr(x)=strcat( reaction.(organ).rxns.abbr(x));
                                % default reaction
                                substr = [param.PeripheralOrgan{i},'_',Metabolite,'[c]'];
                                prod= [param.PeripheralOrgan{i},'_',Metabolite,'[e]'];
                                Organ_secretion_gene = {''}; % should become an input questions
                                %   [model] = addReactionsHH(model, {[param.PeripheralOrgan{i},'_',upper(Metabolite),'t[bc]']},strcat('Transport of ',space,MetaboliteName,' from Cytosol to Extracellular Space (blood) (',param.PeripheralOrgan{i},')'), {[substr,' <==> ',prod]}, {Organ_secretion_gene}, {'Transport, extracellular'},couplingFactor);
                                [model] = addReactionsHH(model, {[param.PeripheralOrgan{i},'_',upper(Metabolite),'t']},strcat('Transport of ',space,MetaboliteName,' from Cytosol to Extracellular Space (blood) (',param.PeripheralOrgan{i},')'), {[substr,' <==> ',prod]}, {Organ_secretion_gene}, {'Transport, extracellular'},couplingFactor);
                            end
                        end
                    end
                    %         model = addOrganRxn(model,reaction, organ,rxnT, couplingFactor,'[bc]');
                    model = addOrganRxn(model,reaction, organ,rxnT, couplingFactor);
                else
                    % default reaction
                    substr = [param.PeripheralOrgan{i},'_',Metabolite,'[c]'];
                    prod= [param.PeripheralOrgan{i},'_',Metabolite,'[e]'];
                    Organ_secretion_gene = {''}; % should become an input questions
                    %  [model] = addReactionsHH(model, {[param.PeripheralOrgan{i},'_',upper(Metabolite),'t[bc]']},strcat('Transport of ',space,MetaboliteName,' from Cytosol to Extracellular Space (blood) (',param.PeripheralOrgan{i},')'), {[substr,' <==> ',prod]}, {Organ_secretion_gene}, {'Transport, extracellular'},couplingFactor);
                    [model] = addReactionsHH(model, {[param.PeripheralOrgan{i},'_',upper(Metabolite),'t']},strcat('Transport of ',space,MetaboliteName,' from Cytosol to Extracellular Space (blood) (',param.PeripheralOrgan{i},')'), {[substr,' <==> ',prod]}, {Organ_secretion_gene}, {'Transport, extracellular'},couplingFactor);
                end
            else
                % default reaction
                substr = [param.PeripheralOrgan{i},'_',Metabolite,'[c]'];
                prod= [param.PeripheralOrgan{i},'_',Metabolite,'[e]'];
                Organ_secretion_gene = {''}; % should become an input questions
                %        [model] = addReactionsHH(model, {[param.PeripheralOrgan{i},'_',upper(Metabolite),'t[bc]']},strcat('Transport of ',space,MetaboliteName,' from Cytosol to Extracellular Space (blood) (',param.PeripheralOrgan{i},')'), {[substr,' <==> ',prod]}, {Organ_secretion_gene}, {'Transport, extracellular'},couplingFactor);
                [model] = addReactionsHH(model, {[param.PeripheralOrgan{i},'_',upper(Metabolite),'t']},strcat('Transport of ',space,MetaboliteName,' from Cytosol to Extracellular Space (blood) (',param.PeripheralOrgan{i},')'), {[substr,' <==> ',prod]}, {Organ_secretion_gene}, {'Transport, extracellular'},couplingFactor);
            end
            % 'Skin_EX_glc_D(e)_[bc]'	'Skin_EX_glc_D(e)_[bc]'	'Skin_glc_D[e]  -> glc_D[bc] ' 'Transport, biofluid'
            substr = [param.PeripheralOrgan{i},'_',Metabolite,'[e]'];
            prod = [Metabolite,'[bc]'];
            [model] = addReactionsHH(model, {[param.PeripheralOrgan{i},'_EX_',Metabolite,'(e)_[bc]']},strcat('Exchange of',space,MetaboliteName,' (',param.PeripheralOrgan{i},')'), {[substr,' <==> ',prod]}, {''}, {'Transport, biofluid'},couplingFactor);
            % add sweat transport
            if isfield(param,'Sweat')  && ~isempty(param.Sweat) && contains(lower(param.Sweat),'y')
                if exist('reaction','var') && isfield(reaction,param.PeripheralOrgan{i})
                    organ = param.PeripheralOrgan{i};
                    % find all relevant reactions
                    e = reaction.(organ).rxns.abbr(contains(reaction.(organ).rxns.formula,'[swS]'));
                    c = reaction.(organ).rxns.abbr(contains(reaction.(organ).rxns.formula,'[c]'));
                    rxnT = intersect(e,c);
                    if ~isempty(rxnT)
                        model = addOrganRxn(model,reaction, organ,rxnT, couplingFactor,'[swS]');
                    else
                        % default reaction
                        substr = [param.PeripheralOrgan{i},'_',Metabolite,'[c]'];
                        prod= [param.PeripheralOrgan{i},'_',Metabolite,'[swS]'];
                        Organ_secretion_gene = {''}; % should become an input questions
                        [model] = addReactionsHH(model, {[param.PeripheralOrgan{i},'_',upper(Metabolite),'t[swS]']},strcat('Transport of ',space,MetaboliteName,' from Cytosol to Extracellular Space (sweat) (',param.PeripheralOrgan{i},')'), {[substr,' --> ',prod]}, {Organ_secretion_gene}, {'Transport, extracellular'},couplingFactor);
                        
                    end
                else
                    % default reaction
                    substr = [param.PeripheralOrgan{i},'_',Metabolite,'[c]'];
                    prod= [param.PeripheralOrgan{i},'_',Metabolite,'[swS]'];
                    Organ_secretion_gene = {''}; % should become an input questions
                    [model] = addReactionsHH(model, {[param.PeripheralOrgan{i},'_',upper(Metabolite),'t[swS]']},strcat('Transport of ',space,MetaboliteName,' from Cytosol to Extracellular Space (sweat) (',param.PeripheralOrgan{i},')'), {[substr,' --> ',prod]}, {Organ_secretion_gene}, {'Transport, extracellular'},couplingFactor);
                end
                %      'Skin_EX_glc_D(swS)_[sw]'
                substr = [param.PeripheralOrgan{i},'_',Metabolite,'[swS]'];
                prod = [Metabolite,'[swS]'];
                [model] = addReactionsHH(model, {[param.PeripheralOrgan{i},'_EX_',Metabolite,'(swS)_[sw]']},strcat('Exchange of',space,MetaboliteName,' (',param.PeripheralOrgan{i},')'), {[substr,' --> ',prod]}, {''}, {'Transport, biofluid'},couplingFactor);
                
                % 'EX_glc_D[sw]'
                substr = [Metabolite,'[sw]'];
                prod = [''];
                [model] = addReactionsHH(model, {['EX_',Metabolite,'[sw]']},strcat('Excretion of',space,MetaboliteName,' (sweat)'), {[substr,' --> ',prod]}, {''}, {'Exchange/demand reaction'},couplingFactor);
            end
        elseif strcmp(param.PeripheralOrgan{i},'Breast')
            
            if exist('reaction','var') && isfield(reaction,param.PeripheralOrgan{i})
                organ = param.PeripheralOrgan{i};
                % transport into the sIEC
                % find all relevant reactions
                e = reaction.(organ).rxns.abbr(contains(reaction.(organ).rxns.formula,'[e]'));
                c = reaction.(organ).rxns.abbr(contains(reaction.(organ).rxns.formula,'[c]'));
                rxnT = intersect(e,c);
                if ~isempty(rxnT)
                    % if this is an abc transporter - also enable import
                    % with default reaction
                    if (contains(rxnT,'ABCt'))
                        % check that indeed an abc efflux
                        for r = 1: length(rxnT)
                            x = find(ismember(reaction.(organ).rxns.abbr,rxnT{r}));
                            if ~isempty(contains(   reaction.(organ).rxns.formula(x),'atp[c]'))
                                % add to abbr [bc]
                                %   reaction.(organ).rxns.abbr(x)=strcat( reaction.(organ).rxns.abbr(x),'[bc]');
                                reaction.(organ).rxns.abbr(x)=strcat( reaction.(organ).rxns.abbr(x));
                                % default reaction
                                substr = [param.PeripheralOrgan{i},'_',Metabolite,'[c]'];
                                prod= [param.PeripheralOrgan{i},'_',Metabolite,'[e]'];
                                Organ_secretion_gene = {''}; % should become an input questions
                                %                  [model] = addReactionsHH(model, {[param.PeripheralOrgan{i},'_',upper(Metabolite),'t[bc]']},strcat('Transport of ',space,MetaboliteName,' from Cytosol to Extracellular Space (blood) (',param.PeripheralOrgan{i},')'), {[substr,' <==> ',prod]}, {Organ_secretion_gene}, {'Transport, extracellular'},couplingFactor);
                                [model] = addReactionsHH(model, {[param.PeripheralOrgan{i},'_',upper(Metabolite),'t']},strcat('Transport of ',space,MetaboliteName,' from Cytosol to Extracellular Space (blood) (',param.PeripheralOrgan{i},')'), {[substr,' <==> ',prod]}, {Organ_secretion_gene}, {'Transport, extracellular'},couplingFactor);
                            end
                        end
                    end
                    %         model = addOrganRxn(model,reaction, organ,rxnT, couplingFactor,'[bc]');
                    model = addOrganRxn(model,reaction, organ,rxnT, couplingFactor);
                else
                    % default reaction
                    substr = [param.PeripheralOrgan{i},'_',Metabolite,'[c]'];
                    prod= [param.PeripheralOrgan{i},'_',Metabolite,'[e]'];
                    Organ_secretion_gene = {''}; % should become an input questions
                    %             [model] = addReactionsHH(model, {[param.PeripheralOrgan{i},'_',upper(Metabolite),'t[bc]']},strcat('Transport of ',space,MetaboliteName,' from Cytosol to Extracellular Space (blood) (',param.PeripheralOrgan{i},')'), {[substr,' <==> ',prod]}, {Organ_secretion_gene}, {'Transport, extracellular'},couplingFactor);
                    [model] = addReactionsHH(model, {[param.PeripheralOrgan{i},'_',upper(Metabolite),'t']},strcat('Transport of ',space,MetaboliteName,' from Cytosol to Extracellular Space (blood) (',param.PeripheralOrgan{i},')'), {[substr,' <==> ',prod]}, {Organ_secretion_gene}, {'Transport, extracellular'},couplingFactor);
                end
            else
                % default reaction
                substr = [param.PeripheralOrgan{i},'_',Metabolite,'[c]'];
                prod= [param.PeripheralOrgan{i},'_',Metabolite,'[e]'];
                Organ_secretion_gene = {''}; % should become an input questions
                %     [model] = addReactionsHH(model, {[param.PeripheralOrgan{i},'_',upper(Metabolite),'t[bc]']},strcat('Transport of ',space,MetaboliteName,' from Cytosol to Extracellular Space (blood) (',param.PeripheralOrgan{i},')'), {[substr,' <==> ',prod]}, {Organ_secretion_gene}, {'Transport, extracellular'},couplingFactor);
                [model] = addReactionsHH(model, {[param.PeripheralOrgan{i},'_',upper(Metabolite),'t']},strcat('Transport of ',space,MetaboliteName,' from Cytosol to Extracellular Space (blood) (',param.PeripheralOrgan{i},')'), {[substr,' <==> ',prod]}, {Organ_secretion_gene}, {'Transport, extracellular'},couplingFactor);
            end
            % 'Breast_EX_glc_D(e)_[bc]'	'Breast_EX_glc_D(e)_[bc]'	'Breast_glc_D[e]  -> glc_D[bc] ' 'Transport, biofluid'
            substr = [param.PeripheralOrgan{i},'_',Metabolite,'[e]'];
            prod = [Metabolite,'[bc]'];
            [model] = addReactionsHH(model, {[param.PeripheralOrgan{i},'_EX_',Metabolite,'(e)_[bc]']},strcat('Exchange of',space,MetaboliteName,' (',param.PeripheralOrgan{i},')'), {[substr,' <==> ',prod]}, {''}, {'Transport, biofluid'},couplingFactor);
            
            if  isfield(param,'BreastMilk')  && ~isempty(param.BreastMilk) && contains(lower(param.BreastMilk),'y')
                if exist('reaction','var') && isfield(reaction,param.PeripheralOrgan{i})
                    organ = param.PeripheralOrgan{i};
                    % transport into the sIEC
                    % find all relevant reactions
                    e = reaction.(organ).rxns.abbr(contains(reaction.(organ).rxns.formula,'[miB]'));
                    c = reaction.(organ).rxns.abbr(contains(reaction.(organ).rxns.formula,'[c]'));
                    rxnT = intersect(e,c);
                    if ~isempty(rxnT)
                        model = addOrganRxn(model,reaction, organ,rxnT, couplingFactor,'[miB]');
                    else
                        % default reaction
                        substr = [param.PeripheralOrgan{i},'_',Metabolite,'[c]'];
                        prod= [param.PeripheralOrgan{i},'_',Metabolite,'[miB]'];
                        Organ_secretion_gene = {''}; % should become an input questions
                        [model] = addReactionsHH(model, {[param.PeripheralOrgan{i},'_',upper(Metabolite),'t[miB]']},strcat('Transport of ',space,MetaboliteName,' from Cytosol to Extracellular Space (breast milk) (',param.PeripheralOrgan{i},')'), {[substr,' --> ',prod]}, {Organ_secretion_gene}, {'Transport, extracellular'},couplingFactor);
                        
                    end
                else
                    % default reaction
                    substr = [param.PeripheralOrgan{i},'_',Metabolite,'[c]'];
                    prod= [param.PeripheralOrgan{i},'_',Metabolite,'[miB]'];
                    Organ_secretion_gene = {''}; % should become an input questions
                    [model] = addReactionsHH(model, {[param.PeripheralOrgan{i},'_',upper(Metabolite),'t[miB]']},strcat('Transport of ',space,MetaboliteName,' from Cytosol to Extracellular Space (breast milk) (',param.PeripheralOrgan{i},')'), {[substr,' --> ',prod]}, {Organ_secretion_gene}, {'Transport, extracellular'},couplingFactor);
                end
                %     Breast_EX_glc_D(miB)_[mi]'
                substr = [param.PeripheralOrgan{i},'_',Metabolite,'[miB]'];
                prod = [Metabolite,'[mi]'];
                [model] = addReactionsHH(model, {[param.PeripheralOrgan{i},'_EX_',Metabolite,'(miB)_[mi]']},strcat('Exchange of',space,MetaboliteName,' (',param.PeripheralOrgan{i},')'), {[substr,' --> ',prod]}, {''}, {'Transport, biofluid'},couplingFactor);
                
                % 'EX_glc_D[mi]'
                substr = [Metabolite,'[mi]'];
                prod = [''];
                [model] = addReactionsHH(model, {['EX_',Metabolite,'[mi]']},strcat('Excretion of',space,MetaboliteName,' (breast milk)'), {[substr,' --> ',prod]}, {''}, {'Exchange/demand reaction'},couplingFactor);
            end
        elseif strcmp(param.PeripheralOrgan{i},'sIEC')
            % sIEC can only take up but not secrete into [bc]
            if exist('reaction','var') && isfield(reaction,param.PeripheralOrgan{i})
                organ = param.PeripheralOrgan{i};
                % transport into the sIEC
                % find all relevant reactions
                e = reaction.(organ).rxns.abbr(contains(reaction.(organ).rxns.formula,'[e]'));
                c = reaction.(organ).rxns.abbr(contains(reaction.(organ).rxns.formula,'[c]'));
                rxnT = intersect(e,c);
                if ~isempty(rxnT)
                    
                    %    model = addOrganRxn(model,reaction, organ,rxnT, couplingFactor,'[bc]');
                    model = addOrganRxn(model,reaction, organ,rxnT, couplingFactor);
                else
                    % default reaction
                    substr = [param.PeripheralOrgan{i},'_',Metabolite,'[c]'];
                    prod= [param.PeripheralOrgan{i},'_',Metabolite,'[e]'];
                    
                    Organ_secretion_gene = {''}; % should become an input questions
                    %                [model] = addReactionsHH(model, {[param.PeripheralOrgan{i},'_',upper(Metabolite),'t[bc]']},strcat('Transport of ',space,MetaboliteName,' from Cytosol to Extracellular Space (blood) (',param.PeripheralOrgan{i},')'), {[substr,' <==> ',prod]}, {Organ_secretion_gene}, {'Transport, extracellular'},couplingFactor);
                    [model] = addReactionsHH(model, {[param.PeripheralOrgan{i},'_',upper(Metabolite),'t']},strcat('Transport of ',space,MetaboliteName,' from Cytosol to Extracellular Space (blood) (',param.PeripheralOrgan{i},')'), {[substr,' <==> ',prod]}, {Organ_secretion_gene}, {'Transport, extracellular'},couplingFactor);
                end
            else % default reaction
                substr = [param.PeripheralOrgan{i},'_',Metabolite,'[c]'];
                prod= [param.PeripheralOrgan{i},'_',Metabolite,'[e]'];
                
                Organ_secretion_gene = {''}; % should become an input questions
                %     [model] = addReactionsHH(model, {[param.PeripheralOrgan{i},'_',upper(Metabolite),'t[bc]']},strcat('Transport of ',space,MetaboliteName,' from Cytosol to Extracellular Space (blood) (',param.PeripheralOrgan{i},')'), {[substr,' <==> ',prod]}, {Organ_secretion_gene}, {'Transport, extracellular'},couplingFactor);
                [model] = addReactionsHH(model, {[param.PeripheralOrgan{i},'_',upper(Metabolite),'t']},strcat('Transport of ',space,MetaboliteName,' from Cytosol to Extracellular Space (blood) (',param.PeripheralOrgan{i},')'), {[substr,' <==> ',prod]}, {Organ_secretion_gene}, {'Transport, extracellular'},couplingFactor);
            end
            % 'Liver_EX_glc_D(e)_[bc]'	'Liver_EX_glc_D(e)_[bc]'	'Liver_glc_D[e]  -> glc_D[bc] ' 'Transport, biofluid'
            substr = [param.PeripheralOrgan{i},'_',Metabolite,'[e]'];
            prod = [Metabolite,'[bc]'];
            [model] = addReactionsHH(model, {[param.PeripheralOrgan{i},'_EX_',Metabolite,'(e)_[bc]']},strcat('Exchange of',space,MetaboliteName,' (',param.PeripheralOrgan{i},')'), {[substr,' <==> ',prod]}, {''}, {'Transport, biofluid'},couplingFactor);
            model.ub(ismember(model.rxns,[param.PeripheralOrgan{i},'_EX_',Metabolite,'(e)_[bc]'])) = 0; % uptake only
            % secretion is into [bp]
            if exist('reaction','var') && isfield(reaction,'sIEC')
                % find all relevant reactions
                e = reaction.(organ).rxns.abbr(contains(reaction.(organ).rxns.formula,'[bpI]'));
                c = reaction.(organ).rxns.abbr(contains(reaction.(organ).rxns.formula,'[c]'));
                rxnT = intersect(e,c);
                if ~isempty(rxnT)
                    model = addOrganRxn(model,reaction, organ,rxnT, couplingFactor,'[bpI]');
                else
                    % default reaction
                    substr = ['sIEC_',Metabolite,'[c]'];
                    prod = ['sIEC_',Metabolite,'[bpI]'];
                    SIEC_secretion_gene = {''}; % should become an input questions
                    [model] = addReactionsHH(model, {['sIEC_',upper(Metabolite),'t[bpI]']},strcat('Transport of ',space,MetaboliteName,' into the portal blood (bpI)'), {[substr,' --> ',prod]}, {SIEC_secretion_gene}, {'Transport, extracellular'},couplingFactor);
                end
            else
                % default reaction
                substr = ['sIEC_',Metabolite,'[c]'];
                prod = ['sIEC_',Metabolite,'[bpI]'];
                SIEC_secretion_gene = {''}; % should become an input questions
                [model] = addReactionsHH(model, {['sIEC_',upper(Metabolite),'t[bpI]']},strcat('Transport of ',space,MetaboliteName,' into the portal blood (bpI)'), {[substr,' --> ',prod]}, {SIEC_secretion_gene}, {'Transport, extracellular'},couplingFactor);
            end
            % excretion into bp
            % 'sIEC_EX_glc_D[bpI]_[bp]'	'Exchange of D-Glucose (sIEC)'	'sIEC_glc_D[bpI]  -> glc_D[bp] '	'Transport, biofluid'
            substr = ['sIEC_',Metabolite,'[bpI]'];
            prod = [Metabolite,'[bp]'];
            [model] = addReactionsHH(model, {['sIEC_EX_',Metabolite,'[bpI]_[bp]']},strcat('Transport of ',space,MetaboliteName,' from Extracellular Space to portal vein (sIEC)'), {[substr,' --> ',prod]}, {''}, {'Transport, biofluid'},couplingFactor);
        elseif strcmp(param.PeripheralOrgan{i},'Colon')
            % Colon can only take up but not secrete into [bc]
            if exist('reaction','var') && isfield(reaction,param.PeripheralOrgan{i})
                organ = param.PeripheralOrgan{i};
                % transport into the Colon
                % find all relevant reactions
                e = reaction.(organ).rxns.abbr(contains(reaction.(organ).rxns.formula,'[e]'));
                c = reaction.(organ).rxns.abbr(contains(reaction.(organ).rxns.formula,'[c]'));
                rxnT = intersect(e,c);
                if ~isempty(rxnT)
                    % if this is an abc transporter - also enable import
                    % with default reaction
                    if (contains(rxnT,'ABCt'))
                        % check that indeed an abc efflux
                        for r = 1: length(rxnT)
                            x = find(ismember(reaction.(organ).rxns.abbr,rxnT{r}));
                            if ~isempty(contains(   reaction.(organ).rxns.formula(x),'atp[c]'))
                                % add to abbr [bc]
                                %           reaction.(organ).rxns.abbr(x)=strcat( reaction.(organ).rxns.abbr(x),'[bc]');
                                reaction.(organ).rxns.abbr(x)=strcat( reaction.(organ).rxns.abbr(x));
                                % default reaction
                                substr = [param.PeripheralOrgan{i},'_',Metabolite,'[c]'];
                                prod= [param.PeripheralOrgan{i},'_',Metabolite,'[e]'];
                                Organ_secretion_gene = {''}; % should become an input questions
                                %              [model] = addReactionsHH(model, {[param.PeripheralOrgan{i},'_',upper(Metabolite),'t[bc]']},strcat('Transport of ',space,MetaboliteName,' from Cytosol to Extracellular Space (blood) (',param.PeripheralOrgan{i},')'), {[substr,' <==> ',prod]}, {Organ_secretion_gene}, {'Transport, extracellular'},couplingFactor);
                                [model] = addReactionsHH(model, {[param.PeripheralOrgan{i},'_',upper(Metabolite),'t']},strcat('Transport of ',space,MetaboliteName,' from Cytosol to Extracellular Space (blood) (',param.PeripheralOrgan{i},')'), {[substr,' <==> ',prod]}, {Organ_secretion_gene}, {'Transport, extracellular'},couplingFactor);
                            end
                        end
                    end
                    %   model = addOrganRxn(model,reaction, organ,rxnT, couplingFactor,'[bc]');
                    model = addOrganRxn(model,reaction, organ,rxnT, couplingFactor);
                else
                    % default reaction
                    substr = [param.PeripheralOrgan{i},'_',Metabolite,'[c]'];
                    prod= [param.PeripheralOrgan{i},'_',Metabolite,'[e]'];
                    
                    Organ_secretion_gene = {''}; % should become an input questions
                    %   [model] = addReactionsHH(model, {[param.PeripheralOrgan{i},'_',upper(Metabolite),'t[bc]']},strcat('Transport of ',space,MetaboliteName,' from Cytosol to Extracellular Space (blood) (',param.PeripheralOrgan{i},')'), {[substr,' <==> ',prod]}, {Organ_secretion_gene}, {'Transport, extracellular'},couplingFactor);
                    [model] = addReactionsHH(model, {[param.PeripheralOrgan{i},'_',upper(Metabolite),'t']},strcat('Transport of ',space,MetaboliteName,' from Cytosol to Extracellular Space (blood) (',param.PeripheralOrgan{i},')'), {[substr,' <==> ',prod]}, {Organ_secretion_gene}, {'Transport, extracellular'},couplingFactor);
                end
            else % default reaction
                substr = [param.PeripheralOrgan{i},'_',Metabolite,'[c]'];
                prod= [param.PeripheralOrgan{i},'_',Metabolite,'[e]'];
                
                Organ_secretion_gene = {''}; % should become an input questions
                %     [model] = addReactionsHH(model, {[param.PeripheralOrgan{i},'_',upper(Metabolite),'t[bc]']},strcat('Transport of ',space,MetaboliteName,' from Cytosol to Extracellular Space (blood) (',param.PeripheralOrgan{i},')'), {[substr,' <==> ',prod]}, {Organ_secretion_gene}, {'Transport, extracellular'},couplingFactor);
                [model] = addReactionsHH(model, {[param.PeripheralOrgan{i},'_',upper(Metabolite),'t']},strcat('Transport of ',space,MetaboliteName,' from Cytosol to Extracellular Space (blood) (',param.PeripheralOrgan{i},')'), {[substr,' <==> ',prod]}, {Organ_secretion_gene}, {'Transport, extracellular'},couplingFactor);
            end
            % 'Liver_EX_glc_D(e)_[bc]'	'Liver_EX_glc_D(e)_[bc]'	'Liver_glc_D[e]  -> glc_D[bc] ' 'Transport, biofluid'
            substr = [param.PeripheralOrgan{i},'_',Metabolite,'[e]'];
            prod = [Metabolite,'[bc]'];
            [model] = addReactionsHH(model, {[param.PeripheralOrgan{i},'_EX_',Metabolite,'(e)_[bc]']},strcat('Exchange of',space,MetaboliteName,' (',param.PeripheralOrgan{i},')'), {[substr,' <==> ',prod]}, {''}, {'Transport, biofluid'},couplingFactor);
            model.ub(ismember(model.rxns,[param.PeripheralOrgan{i},'_EX_',Metabolite,'(e)_[bc]'])) = 0; % uptake only
            % secretion is into [bp]
            if exist('reaction','var') && isfield(reaction,'Colon')
                % find all relevant reactions
                e = reaction.(organ).rxns.abbr(contains(reaction.(organ).rxns.formula,'[bpC]'));
                c = reaction.(organ).rxns.abbr(contains(reaction.(organ).rxns.formula,'[c]'));
                rxnT = intersect(e,c);
                if ~isempty(rxnT)
                    reactionTmp = reaction;
                    reactionTmp.(organ).rxns.formula = regexprep(reaction.(organ).rxns.formula,'\[e\]','[bpC]');
                    model = addOrganRxn(model,reaction, organ,rxnT, couplingFactor,'[bpC]');
                    clear reactionTmp
                else
                    % default reaction
                    substr = ['Colon_',Metabolite,'[c]'];
                    prod = ['Colon_',Metabolite,'[bpC]'];
                    SIEC_secretion_gene = {''}; % should become an input questions
                    [model] = addReactionsHH(model, {['Colon_',upper(Metabolite),'t[bpC]']},strcat('Transport of ',space,MetaboliteName,' into the portal blood (bpC)'), {[substr,' --> ',prod]}, {SIEC_secretion_gene}, {'Transport, extracellular'},couplingFactor);
                end
            else
                % default reaction
                substr = ['Colon_',Metabolite,'[c]'];
                prod = ['Colon_',Metabolite,'[bpC]'];
                SIEC_secretion_gene = {''}; % should become an input questions
                [model] = addReactionsHH(model, {['Colon_',upper(Metabolite),'t[bpC]']},strcat('Transport of ',space,MetaboliteName,' into the portal blood (bpC)'), {[substr,' --> ',prod]}, {SIEC_secretion_gene}, {'Transport, extracellular'},couplingFactor);
            end
            % excretion into bp
            % 'sIEC_EX_glc_D[bpI]_[bp]'	'Exchange of D-Glucose (sIEC)'	'sIEC_glc_D[bpI]  -> glc_D[bp] '	'Transport, biofluid'
            substr = ['Colon_',Metabolite,'[bpC]'];
            prod = [Metabolite,'[bp]'];
            [model] = addReactionsHH(model, {['Colon_EX_',Metabolite,'[bpC]_[bp]']},strcat('Transport of ',space,MetaboliteName,' from Extracellular Space to portal vein (Colon)'), {[substr,' --> ',prod]}, {''}, {'Transport, biofluid'},couplingFactor);
            
        elseif strcmp(param.PeripheralOrgan{i},'Stomach')
            %% uptake by stomach (currently only bc exchange)
            if exist('reaction','var') && isfield(reaction,param.PeripheralOrgan{i})
                organ = param.PeripheralOrgan{i};
                % transport into the sIEC
                % find all relevant reactions
                e = reaction.(organ).rxns.abbr(contains(reaction.(organ).rxns.formula,'[e]'));
                c = reaction.(organ).rxns.abbr(contains(reaction.(organ).rxns.formula,'[c]'));
                rxnT = intersect(e,c);
                if ~isempty(rxnT)
                    % if this is an abc transporter - also enable import
                    % with default reaction
                    if (contains(rxnT,'ABCt'))
                        % check that indeed an abc efflux
                        for r = 1: length(rxnT)
                            x = find(ismember(reaction.(organ).rxns.abbr,rxnT{r}));
                            if ~isempty(contains(   reaction.(organ).rxns.formula(x),'atp[c]'))
                                % add to abbr [bc]
                                %  reaction.(organ).rxns.abbr(x)=strcat( reaction.(organ).rxns.abbr(x),'[bc]');
                                reaction.(organ).rxns.abbr(x)=strcat( reaction.(organ).rxns.abbr(x));
                                % default reaction
                                substr = [param.PeripheralOrgan{i},'_',Metabolite,'[c]'];
                                prod= [param.PeripheralOrgan{i},'_',Metabolite,'[e]'];
                                Organ_secretion_gene = {''}; % should become an input questions
                                %       [model] = addReactionsHH(model, {[param.PeripheralOrgan{i},'_',upper(Metabolite),'t[bc]']},strcat('Transport of ',space,MetaboliteName,' from Cytosol to Extracellular Space (blood) (',param.PeripheralOrgan{i},')'), {[substr,' <==> ',prod]}, {Organ_secretion_gene}, {'Transport, extracellular'},couplingFactor);
                                [model] = addReactionsHH(model, {[param.PeripheralOrgan{i},'_',upper(Metabolite),'t']},strcat('Transport of ',space,MetaboliteName,' from Cytosol to Extracellular Space (blood) (',param.PeripheralOrgan{i},')'), {[substr,' <==> ',prod]}, {Organ_secretion_gene}, {'Transport, extracellular'},couplingFactor);
                            end
                        end
                    end
                    %           model = addOrganRxn(model,reaction, organ,rxnT, couplingFactor,'[bc]');
                    model = addOrganRxn(model,reaction, organ,rxnT, couplingFactor);
                else
                    % default reaction
                    substr = [param.PeripheralOrgan{i},'_',Metabolite,'[c]'];
                    prod= [param.PeripheralOrgan{i},'_',Metabolite,'[e]'];
                    Organ_secretion_gene = {''}; % should become an input questions
                    %    [model] = addReactionsHH(model, {[param.PeripheralOrgan{i},'_',upper(Metabolite),'t[bc]']},strcat('Transport of ',space,MetaboliteName,' from Cytosol to Extracellular Space (blood) (',param.PeripheralOrgan{i},')'), {[substr,' <=> ',prod]}, {Organ_secretion_gene}, {'Transport, extracellular'},couplingFactor);
                    [model] = addReactionsHH(model, {[param.PeripheralOrgan{i},'_',upper(Metabolite),'t']},strcat('Transport of ',space,MetaboliteName,' from Cytosol to Extracellular Space (blood) (',param.PeripheralOrgan{i},')'), {[substr,' <=> ',prod]}, {Organ_secretion_gene}, {'Transport, extracellular'},couplingFactor);
                end
            else
                % default reaction
                substr = [param.PeripheralOrgan{i},'_',Metabolite,'[c]'];
                prod= [param.PeripheralOrgan{i},'_',Metabolite,'[e]'];
                Organ_secretion_gene = {''}; % should become an input questions
                %    [model] = addReactionsHH(model, {[param.PeripheralOrgan{i},'_',upper(Metabolite),'t[bc]']},strcat('Transport of ',space,MetaboliteName,' from Cytosol to Extracellular Space (blood) (',param.PeripheralOrgan{i},')'), {[substr,' <==> ',prod]}, {Organ_secretion_gene}, {'Transport, extracellular'},couplingFactor);
                [model] = addReactionsHH(model, {[param.PeripheralOrgan{i},'_',upper(Metabolite),'t']},strcat('Transport of ',space,MetaboliteName,' from Cytosol to Extracellular Space (blood) (',param.PeripheralOrgan{i},')'), {[substr,' <==> ',prod]}, {Organ_secretion_gene}, {'Transport, extracellular'},couplingFactor);
            end
            % 'Breast_EX_glc_D(e)_[bc]'	'Breast_EX_glc_D(e)_[bc]'	'Breast_glc_D[e]  -> glc_D[bc] ' 'Transport, biofluid'
            substr = [param.PeripheralOrgan{i},'_',Metabolite,'[e]'];
            prod = [Metabolite,'[bc]'];
            [model] = addReactionsHH(model, {[param.PeripheralOrgan{i},'_EX_',Metabolite,'(e)_[bc]']},strcat('Exchange of',space,MetaboliteName,' (',param.PeripheralOrgan{i},')'), {[substr,' <==> ',prod]}, {''}, {'Transport, biofluid'},couplingFactor);
            
            if  isfield(param,'StomachUptake')  && ~isempty(param.StomachUptake) && contains(lower(param.StomachUptake),'y')
                
                %     Breast_EX_glc_D(miB)_[mi]'
                substr = [param.PeripheralOrgan{i},'_',Metabolite,'[luS]'];
                prod = [Metabolite,'[lu]'];
                [model] = addReactionsHH(model, {[param.PeripheralOrgan{i},'_EX_',Metabolite,'(luS)_[lu]']},strcat('Exchange of',space,MetaboliteName,' (',param.PeripheralOrgan{i},')'), {[substr,' --> ',prod]}, {''}, {'Transport, biofluid'},couplingFactor);
                
            end
        else
            
            if exist('reaction','var') && isfield(reaction,param.PeripheralOrgan{i})
                organ = param.PeripheralOrgan{i};
                % transport into the organ
                e = reaction.(organ).rxns.abbr(contains(reaction.(organ).rxns.formula,'[e]'));
                c = reaction.(organ).rxns.abbr(contains(reaction.(organ).rxns.formula,'[c]'));
                rxnT = intersect(e,c);
                
                if ~isempty(rxnT)
                    % this does not mean that there is a transport reaction for
                    % the metabolite in question
                    x = find(ismember(reaction.(organ).rxns.abbr,rxnT));
                    if contains(reaction.(organ).rxns.formula(x),strcat(Metabolite,'['))
                        % if this is an abc transporter - also enable import
                        % with default reaction
                        if (contains(rxnT,'ABCt'))
                            % check that indeed an abc efflux
                            for r = 1: length(rxnT)
                                x = find(ismember(reaction.(organ).rxns.abbr,rxnT{r}));
                                if ~isempty(contains(   reaction.(organ).rxns.formula(x),'atp[c]'))
                                    % add to abbr [bc]
                                    % reaction.(organ).rxns.abbr(x)=strcat( reaction.(organ).rxns.abbr(x),'[bc]');
                                    reaction.(organ).rxns.abbr(x)=strcat( reaction.(organ).rxns.abbr(x));
                                    % default reaction
                                    substr = [param.PeripheralOrgan{i},'_',Metabolite,'[c]'];
                                    prod= [param.PeripheralOrgan{i},'_',Metabolite,'[e]'];
                                    Organ_secretion_gene = {''}; % should become an input questions
                                    %   [model] = addReactionsHH(model, {[param.PeripheralOrgan{i},'_',upper(Metabolite),'t[bc]']},strcat('Transport of ',space,MetaboliteName,' from Cytosol to Extracellular Space (blood) (',param.PeripheralOrgan{i},')'), {[substr,' <==> ',prod]}, {Organ_secretion_gene}, {'Transport, extracellular'},couplingFactor);
                                    [model] = addReactionsHH(model, {[param.PeripheralOrgan{i},'_',upper(Metabolite),'t']},strcat('Transport of ',space,MetaboliteName,' from Cytosol to Extracellular Space (blood) (',param.PeripheralOrgan{i},')'), {[substr,' <==> ',prod]}, {Organ_secretion_gene}, {'Transport, extracellular'},couplingFactor);
                                end
                            end
                        end
                        %  model = addOrganRxn(model,reaction, organ,rxnT, couplingFactor,'[bc]');
                        model = addOrganRxn(model,reaction, organ,rxnT, couplingFactor);
                    else
                        % default reaction
                        substr = [param.PeripheralOrgan{i},'_',Metabolite,'[c]'];
                        prod= [param.PeripheralOrgan{i},'_',Metabolite,'[e]'];
                        
                        Organ_secretion_gene = {''}; % should become an input questions
                        %    [model] = addReactionsHH(model, {[param.PeripheralOrgan{i},'_',upper(Metabolite),'t[bc]']},strcat('Transport of ',space,MetaboliteName,' from Cytosol to Extracellular Space (blood) (',param.PeripheralOrgan{i},')'), {[substr,' <==> ',prod]}, {Organ_secretion_gene}, {'Transport, extracellular'},couplingFactor);
                        [model] = addReactionsHH(model, {[param.PeripheralOrgan{i},'_',upper(Metabolite),'t']},strcat('Transport of ',space,MetaboliteName,' from Cytosol to Extracellular Space (blood) (',param.PeripheralOrgan{i},')'), {[substr,' <==> ',prod]}, {Organ_secretion_gene}, {'Transport, extracellular'},couplingFactor);
                        
                    end
                else
                    % default reaction
                    substr = [param.PeripheralOrgan{i},'_',Metabolite,'[c]'];
                    prod= [param.PeripheralOrgan{i},'_',Metabolite,'[e]'];
                    
                    Organ_secretion_gene = {''}; % should become an input questions
                    %         [model] = addReactionsHH(model, {[param.PeripheralOrgan{i},'_',upper(Metabolite),'t[bc]']},strcat('Transport of ',space,MetaboliteName,' from Cytosol to Extracellular Space (blood) (',param.PeripheralOrgan{i},')'), {[substr,' <==> ',prod]}, {Organ_secretion_gene}, {'Transport, extracellular'},couplingFactor);
                    [model] = addReactionsHH(model, {[param.PeripheralOrgan{i},'_',upper(Metabolite),'t']},strcat('Transport of ',space,MetaboliteName,' from Cytosol to Extracellular Space (blood) (',param.PeripheralOrgan{i},')'), {[substr,' <==> ',prod]}, {Organ_secretion_gene}, {'Transport, extracellular'},couplingFactor);
                end
            else % default reaction
                substr = [param.PeripheralOrgan{i},'_',Metabolite,'[c]'];
                prod= [param.PeripheralOrgan{i},'_',Metabolite,'[e]'];
                
                Organ_secretion_gene = {''}; % should become an input questions
                %      [model] = addReactionsHH(model, {[param.PeripheralOrgan{i},'_',upper(Metabolite),'t[bc]']},strcat('Transport of ',space,MetaboliteName,' from Cytosol to Extracellular Space (blood) (',param.PeripheralOrgan{i},')'),{[substr,' <==> ',prod]}, {Organ_secretion_gene}, {'Transport, extracellular'},couplingFactor);
                [model] = addReactionsHH(model, {[param.PeripheralOrgan{i},'_',upper(Metabolite),'t']},strcat('Transport of ',space,MetaboliteName,' from Cytosol to Extracellular Space (blood) (',param.PeripheralOrgan{i},')'),{[substr,' <==> ',prod]}, {Organ_secretion_gene}, {'Transport, extracellular'},couplingFactor);
            end
            % 'Liver_EX_glc_D(e)_[bc]'	'Liver_EX_glc_D(e)_[bc]'	'Liver_glc_D[e]  -> glc_D[bc] ' 'Transport, biofluid'
            substr = [param.PeripheralOrgan{i},'_',Metabolite,'[e]'];
            prod = [Metabolite,'[bc]'];
            [model] = addReactionsHH(model, {[param.PeripheralOrgan{i},'_EX_',Metabolite,'(e)_[bc]']},strcat('Exchange of',space,MetaboliteName,' (',param.PeripheralOrgan{i},')'), {[substr,' <==> ',prod]}, {''}, {'Transport, biofluid'},couplingFactor);
        end
        
    end
    
end
%% add reactions from reaction structure (except for transporters that should have been added beforehand
if exist('reaction','var') && ~isempty(reaction)
    Organs = fieldnames(reaction);
    for i = 1 : length(Organs)
        Organ = Organs{i};
        %         remove = zeros(length(reaction.(Organ).rxns.subSystem),1);
        %         for j = 1 : length(reaction.(Organ).rxns.subSystem)
        %             if find(contains(reaction.(Organ).rxns.subSystem{j},'Transport'))
        %                 remove(j,1) = 1;
        %             end
        %         end
        %         if ~isempty(remove)
        %             remove=find(remove);
        %             reaction.(Organ).rxns.abbr(remove) = [];
        %             reaction.(Organ).rxns.name(remove) = [];
        %             reaction.(Organ).rxns.GPR(remove) = [];
        %             reaction.(Organ).rxns.formula(remove) = [];
        %             reaction.(Organ).rxns.subSystem(remove) = [];
        %         end
        %model.rxns = regexprep(model.rxns,'\[bc\]\[bc\]','[bc]');% fix this issue
        for j = 1 : length(reaction.(Organ).rxns.abbr)
            clear reactionTmp
            if ~isempty(reaction.(Organ).rxns.abbr{j}) %&& ~contains(lower(reaction.(Organ).rxns.subSystem{j}), 'transport')
                i
                j
                reactionTmp.(Organ).rxns.abbr = reaction.(Organ).rxns.abbr(j);
                reactionTmp.(Organ).rxns.name = reaction.(Organ).rxns.name(j);
                reactionTmp.(Organ).rxns.GPR = reaction.(Organ).rxns.GPR(j);
                reactionTmp.(Organ).rxns.formula = reaction.(Organ).rxns.formula(j);
                reactionTmp.(Organ).rxns.subSystem = reaction.(Organ).rxns.subSystem{j};
                model = addOrganRxn(model,reactionTmp, Organ,'', couplingFactor);
                %model.rxns = regexprep(model.rxns,'\[bc\]\[bc\]','[bc]');% fix this issue
            end
        end
    end
end
if isfield(param,'sinks') && ~isempty(param.sinks)
    if isfield(param.sinks,'organs') && ~isempty(param.sinks.organs)
        % adds sink reactions to the listed organs % assumes [c]
        for i = 1 : length(param.sinks.organs)
            met = strcat(param.sinks.organs{i},Metabolite,'[c]');
            model = addSinkReactions(model,met);
            clear met
        end
    end
    if isfield(param.sinks,'biofluids') && ~isempty(param.sinks.biofluids)
        % adds sink reactions to the listed organs % assumes [c]
        for i = 1 : length(param.sinks.biofluids)
            met = strcat(Metabolite,param.sinks.biofluids{i});
            model = addSinkReactions(model,met);
            clear met
        end
    end
end

%model.rxns = regexprep(model.rxns,'\[bc\]\[bc\]','[bc]');% fix this issue
function model = addOrganRxn(model,reaction, Organ, abbr, couplingFactor,biofluid)
if ~exist('biofluid','var')
    biofluid = '';
end
if ~exist('abbr','var') || isempty(abbr)
    rxnAbbr = strcat(reaction.(Organ).rxns.abbr,biofluid);
    rxnName = reaction.(Organ).rxns.name;
    rxnGPR = reaction.(Organ).rxns.GPR;
    rxnForm = reaction.(Organ).rxns.formula;
    rxnSS = reaction.(Organ).rxns.subSystem;
else
    r = find(ismember(reaction.(Organ).rxns.abbr,abbr));
    if isempty(r)
        r = find(ismember(reaction.(Organ).rxns.abbr,strcat(abbr,biofluid)));
        rxnAbbr = (reaction.(Organ).rxns.abbr(r));
    else
        rxnAbbr = strcat(reaction.(Organ).rxns.abbr(r),biofluid);
    end
    rxnName = reaction.(Organ).rxns.name(r);
    rxnGPR = reaction.(Organ).rxns.GPR(r);
    rxnForm = reaction.(Organ).rxns.formula(r);
    rxnSS = reaction.(Organ).rxns.subSystem(r);
end
% if strcmp(biofluid,'[bc]')
%     reaction.(Organ).rxns.abbr = strcat(reaction.(Organ).rxns.abbr,'[bc]');
% % elseif strcmp(biofluid,'[bpC]')
% %     rxnAbbr = strcat(rxnAbbr,'[bpC]');
% end
rxnAbbr = regexprep(rxnAbbr,'\[bc\]\[bc\]','[bc]');% fix this issue
rxnAbbr = regexprep(rxnAbbr,'\[bc\]\[bpC\]','[bc]');% fix this issue
% if ~isempty(biofluid) && ~contains(biofluid,'bc')
%     % replace [e] with biofluid
%          rxnForm = regexprep(rxnForm,'\[e\]',biofluid);
% end
if~isempty(rxnAbbr)
    for x = 1 : length(rxnAbbr)
        if isempty(find(ismember(model.rxns,rxnAbbr{x}),1))
            rxnAbbr(x);
            [model] = addReactionsHH(model,rxnAbbr(x),rxnName(x), rxnForm(x), rxnGPR(x), rxnSS(x),couplingFactor);
        end
    end
end