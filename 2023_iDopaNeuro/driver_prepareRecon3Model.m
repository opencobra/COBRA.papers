% Recon 3D
% Preparation of the generic Recon3D model with thermodynamic data, and the
% addition of names for reactions where it is missing. It excludes the 
% expression of the reactions because that isspecific information which it 
% is included in XomicsToModel.
%
% The generic model is saved in :
% '~/work/sbgCloud/programExperimental/projects/xomics/data/Recon3D_301/Recon3DModel_301_xomics_input.mat'

clear

inputFolder = ['~' filesep 'work' filesep 'sbgCloud' filesep 'programExperimental' filesep 'projects' filesep 'xomics' filesep 'data' filesep 'Recon3D_301'];
genericModelName = 'Recon3DModel_301_xomics_input';

if 0
    %apply vonBertalanffy 2.* to estimate thermochemical parameters
    run ~/work/sbgCloud/code/fork-COBRA.tutorials/analysis/vonBertalanffy/tutorial_vonBertalanffy_Recon3DModel_301.mlx
else
    %copy the result to the xomics data folder
    copyfile('~/work/sbgCloud/programExperimental/projects/xomics/data/Recon3D_301/20210520T123849_Recon3DModel_301_thermo.mat',[inputFolder filesep 'Recon3DModel_301_thermo.mat'])
end



if isfile([inputFolder filesep genericModelName '.mat'])
    display('Generic model already prepared...')
else
    display('Preparing generic model ...')
    
    %Recon3DModel_301_xomics_input is the output of xomics/code/scripts/driver_prepareRecon3Model.m
    load([inputFolder filesep 'Recon3DModel_301_thermo.mat'])
    
    % model.DrGt0 - `n x 1` array of standard transformed reaction Gibbs energies in kJ/mol.
    if ~isfield(model,'DrGt0')
        error([inputFolder filesep 'Recon3DModel_301_thermo.mat is missing thermochemical data'])
    else
        if isfield(model, 'Srecon')
            model.Sthermo = model.S;
            model.S = model.Srecon;
            model = rmfield(model,'Srecon');
        end
        if isfield(model, 'DfG0_pseudoisomers')
            model = rmfield(model,'DfG0_pseudoisomers');
        end
        if isfield(model, 'directions')
            model = rmfield(model,'directions');
        end
        if isfield(model, 'DfG0_cc_cov')
            model = rmfield(model, 'DfG0_cc_cov');
        end
        if isfield(model, 'DrG0_cc_cov')
            model = rmfield(model, 'DrG0_cc_cov');
        end
        modelThermo = model;
    end
    
    % RXN files directory (for bond enthaplies and bonds broken and formed)
    rxnDir = ['~' filesep 'work' filesep 'code' filesep 'ctf' filesep 'rxns' filesep 'atomMapped'];
    if ~isfolder(rxnDir)
        
        %this needs to be fixed so that the complete input can be generated
        load([inputFolder filesep 'Recon3DModel_301_thermo_BBF.mat'])
        modelThermo.bondsBF = model.bondsBF';
        modelThermo.bondsE = model.bondsE';
        if 1
            modelThermo.meanBBF =  mean(modelThermo.bondsBF(modelThermo.SConsistentRxnBool & modelThermo.bondsBF ~=0 ), 'omitnan');
            modelThermo.meanBE =  mean(modelThermo.bondsBF(modelThermo.SConsistentRxnBool & modelThermo.bondsBF ~=0 ) ,'omitnan');
        else
            modelThermo.meanBBF = model.meanBBF;
            modelThermo.meanBE = model.meanBE;
        end
        modelThermo.SConsistentRxnBool = model.SConsistentRxnBool;
        modelThermo.SConsistentMetBool = model.SConsistentMetBool;
        
        model = modelThermo;
        clear modelThermo;
        
    else
        
        %extracted from Recon3D_BE_BBF_RxnExp.m
        if isfield(model, 'bondsBF')
            model = rmfield(model,'bondsBF');
        end
        if isfield(model, 'bondsE')
            model = rmfield(model,'bondsE');
        end
        if isfield(model, 'meanBBF')
            model = rmfield(model,'meanBBF');
        end
        if isfield(model, 'meanBE')
            model = rmfield(model,'meanBE');
        end
        
        %% Get bond enthalpies and bonds broken and formed
        % Calculate bond enthalpies and number of bonds broken and formed
        [bondsBF, bondsE, meanBBF, meanBE] = findBEandBBF(model, rxnDir, 1);
        
        % Add data to the model
        model.bondsBF = bondsBF;
        model.bondsE = bondsE;
        model.meanBBF = meanBBF;
        model.meanBE = meanBE;
        
    end
    
    if 0
        % Gene expression is added in the XomicsToModel function
        % Generate expressionRxns
        dataFolder = ['~' filesep 'work' filesep 'sbgCloud' filesep 'programReconstruction' ...
            filesep 'projects' filesep 'exoMetDN' filesep 'data' filesep 'xomics' filesep];
        transcriptomicFile = 'SM7.tab';
        transcriptomicData = transcriptomics2genes([dataFolder transcriptomicFile]);
        [model.expressionRxns, ~] = mapExpressionToReactions(model, transcriptomicData);
    end
    
    if 0
        % No lumped reactions is less accurate
        remove lumped reactions
        [row,col,s] = find(model.S >2 | model.S <-2);
        lumpedBool=false(size(model.S,2),1);
        lumpedBool(col)=1;
        transportRxnBool = transportReactionBool(model);
        lumpedBool = lumpedBool & ~transportRxnBool;
        model = removeRxns(model, model.rxns(lumpedBool), 'metRemoveMethod', 'exclusive', 'ctrsRemoveMethod','infeasible');
    end
    
    % Replace missing names (moved from XomicsToModel function)
    model.rxnNames{ismember(model.rxns, 'EX_M03117[e]')} = 'Exchange of Undecanoic acid';
    model.rxnNames{ismember(model.rxns, 'EX_icit[e]')} = 'Exchange of Isocitrate';
    %% Save the model to be input into xomics2model
    
    save([inputFolder filesep genericModelName], 'model')
end

display(['Generic model is: ' inputFolder filesep genericModelName '.mat'])

