%Compare the properties of different human metabolic reconstruction versions for the cardinality optimisation paper

%choose the input data
%assumes the input data is already within the directory
%modelCollectionDirectory='~/work/sbgCloud/programReconstruction/projects/recon2models/data/reconXComparisonModels/';
modelCollectionDirectory='~/work/sbgCloud/programModelling/projects/cardinalityOpt/data/humanComparisonModels/';
cd(modelCollectionDirectory)
    
%choose results directory (without the / at the end!!)
resultsDirectory='~/work/sbgCloud/programModelling/projects/cardinalityOpt/results/humanComparison';

modelConsistencyPlot =1;

%choose the printLevel
printLevel=2;
%
resultsBaseFileName=['humanComparison_results_' datestr(now,30)];
resultsFileName=[resultsBaseFileName '.mat'];

if 0
    approach='single';
else
    approach='batch';
end
%approach='resultsTable';

% modelMetaData         Cell array, where each row is metadata for one model
%                       with five columns: species, modelID, fileName, PMID, doi.
%                       See function modelMetaData=modelCitations() for
%                       example. Table columns ordered by order of rows in
%                       modelMetaData.

%model metadata used to organise and annotate the table of model properties
modelMetaData={...
'Homo sapiens','Recon1.0','Recon1','17267599','';...%version
%'Homo sapiens','Recon2.0','Recon2','23455439','';...%version
%'Homo sapiens','Recon2.0model','Recon2model','23455439','';...%version
%'Homo sapiens','Recon2.02','Recon2v02','23455439','';...%version
%'Homo sapiens','Recon2.02model','Recon2v02model','23455439','';...%version
'Homo sapiens','Recon2.04model','Recon2v04model','23455439','';...%version
'Homo sapiens','Recon2.2model','Recon2v2model','27358602','';...%version
'Homo sapiens','HMR2.0','HMRdatabase2','26209309','';...%version
'Homo sapiens','Recon3.01','Recon301','','';...%version
%'Homo sapiens','Recon3D','Recon301','','';...%version
'Homo sapiens','Recon3.01model','Recon301model','','';...%version
'Homo sapiens','Human1.0','humanGEM1','','';...% v1.0.0 https://github.com/SysBioChalmers/Human-GEM/releases/tag/v1.0.0 12 Mar 2019
%'Homo sapiens','humanGEM1','humanGEM1','','';...% v1.0.0 https://github.com/SysBioChalmers/Human-GEM/releases/tag/v1.0.0 12 Mar 2019
%'Homo sapiens','Human1.10','humanGEM1p10','','';...%https://github.com/SysBioChalmers/Human-GEM/tree/main/model sha 26005b7 16 Sept 2021
%'Homo sapiens','humanGEM1p10','humanGEM1p10','','';...%https://github.com/SysBioChalmers/Human-GEM/tree/main/model sha 26005b7 16 Sept 2021
%'Homo sapiens','Harvetta','Harvetta','','';...
%'Homo sapiens','Harvey','Harvey','','';...
};

[solverOK] = changeCobraSolver('gurobi','LP'); 
[solverOK] = changeCobraSolver('gurobi','QP');

%%
if ~exist(resultsDirectory,'file')
    mkdir(resultsDirectory);
end
cd(resultsDirectory)

switch approach
    %by default, it will run all models in modelCollectionDirectory, but here
    %you can set it up to run only one model
    case 'single'
        %load single model
        if 1
            modelID='Recon1.0';
            load([modelCollectionDirectory modelID '.mat'])
            model=Recon1;
        end
        if 0
            modelID='Recon2.02';
            load([modelCollectionDirectory modelID '.mat'])
            model=Recon202;
        end
        if 0
            modelID='Recon3.01';
            load([modelCollectionDirectory modelID '.mat'])
            model=Recon301;
        end
        if 0
            modelID='humanGEM1';
            load([modelCollectionDirectory modelID '.mat'])
            model=humanGEM1;
        end
        
        if 0
            modelID='HMR2.0';
            load([modelCollectionDirectory modelID '.mat'])
            model=humanGEM1;
        end
        %%%%%%%%%%%
        model = checkModelProperties(model,printLevel);
        %%%%%%%%%%%%%
        k=1;
        modelResults(k).model=model;
        
        [modelResultsTable,modelResults]=makeModelPropertiesTable(modelResults,modelMetaData);
        
        if modelConsistencyPlot
            schematicFlag=0;
            nRows=1;
            nCols=1;
            figureFileName=resultsBaseFileName;
            plotModelConsistency(modelResults,modelMetaData,schematicFlag,nRows,nCols,resultsDirectory,figureFileName);
        end
        
        %results filename timestamped
        save([resultsDirectory filesep resultsFileName],'modelResults','modelMetaData','modelResultsTable');
    case 'batch'
        if 1
            %save results individually for each model
            cd(modelCollectionDirectory);
            %loads the input data into the above directory
            loadHumanModels
        end
        %batch of models in .mat format in a directory
        %assumes that each .mat file is a model
        matFiles=dir(modelCollectionDirectory);
        matFiles.name;
        
        %get rid of entries that do not have a .mat suffix
        bool=false(length(matFiles),1);
        for k=3:length(matFiles)
            if strcmp(matFiles(k).name(end-3:end),'.mat')
                bool(k)=1;
            end
        end
        matFiles=matFiles(bool);
        %     matFiles.name
        
        if 0
            %checks if the models can be loaded and if some exchange reactions
            %can be identified
            for k=1:length(matFiles)
                disp(k)
                disp(matFiles(k).name)
                whosFile=whos('-file',matFiles(k).name);
                if ~strcmp(matFiles(k).name,'clone1.log')
                    load(matFiles(k).name);
                    model=eval(whosFile.name);
                    model=findSExRxnInd(model);
                    printLevel=1;
                end
            end
        end
        

        %generate the results for each mat file, unless the results already exist
        for k=1:length(matFiles)
            resultsFileNamePrefix='modelResults_';
            %create the results structure
            modelResults=struct();
            modelResults(k).matFile=matFiles(k);
            fprintf('%u\t%s\n',k,matFiles(k).name)
            whosFile=whos('-file',matFiles(k).name);
            modelResults(k).modelFilename=matFiles(k).name;
            modelResults(k).modelID=modelResults(k).modelFilename(1:end-4);%take off .mat from filename
            
            if exist([resultsDirectory filesep resultsFileNamePrefix modelResults(k).modelID '.mat'],'file')
                fprintf('%s\n',[resultsDirectory filesep resultsFileNamePrefix modelResults(k).modelID '.mat results found.'])
            else
                fprintf('%s\n',[resultsDirectory filesep resultsFileNamePrefix modelResults(k).modelID '.mat results not found. Being generated...'])
                %load the input data
                load(matFiles(k).name);
                model=eval(whosFile.name);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                model = checkModelProperties(model,printLevel);
                modelResults(k).model=model;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if printLevel>0 && model.FRrowRankDeficiency>0 && 0
                    filePathName=[resultsDirectory  filesep resultsFileNamePrefix modelResults(k).modelID '_rowDependencies.txt'];
                    fileID = fopen(filePathName,'a');
                    fprintf(fileID,'%s\n',matFiles(k).name)
                    fclose(fileID);
                    printFRdependencies(model,filePathName);
                    fileID = fopen(filePathName,'a');
                    fprintf(fileID,'%s\n','-------------------------------')
                    fclose(fileID);
                end
                
                save([resultsDirectory filesep resultsFileNamePrefix modelResults.modelID '.mat'],'modelResults','modelMetaData','-v7.3');
                fprintf('%s\n',[resultsDirectory filesep resultsFileNamePrefix modelResults(k).modelID '.mat results saved.'])
                clear modelResults model;
            end
        end

        %% Amalgamate the results from each model
        %find the names of each modelResults .mat file
        matFiles=dir(resultsDirectory);
        matFiles.name
        %get rid of entries that do not have a .mat suffix
        bool=false(length(matFiles),1);
        for k=3:length(matFiles)
            if strcmp(matFiles(k).name(end-3:end),'.mat') && isempty(strfind(matFiles(k).name,'results'))
                bool(k)=1;
            end
        end
        matFiles=matFiles(bool);
        
        %model properties structure
        modelResultsTmp=struct();
        for k=1:length(matFiles)
            load([resultsDirectory filesep matFiles(k).name])
            tmp=strrep(matFiles(k).name,'modelResults_','');
            fprintf('%u\t%s\n',k,tmp(1:end-4))
            %amalgamate all results into one structure
            modelResultsTmp(k).model=modelResults(k).model;
        end
        clear modelResults;
        modelResults=modelResultsTmp;
        clear modelResultsTmp;
        
        %metadata about each model
        if ~exist('modelMetaData','var')
            if 1
                modelMetaData={'testModel','testModel',modelResults(k).model.modelID,'testModel','testModel'};
            else
                %citations data for models
                modelMetaData=modelCitations();
            end
        end
        
        if 1
            [modelResultsTable,modelResults]=makeModelPropertiesTable(modelResults,modelMetaData);
        end
        
        if modelConsistencyPlot
            schematicFlag=1;
            nRows=4;
            nCols=2;
            figureFileName=resultsBaseFileName;
            plotModelConsistency(modelResults,modelMetaData,schematicFlag,nRows,nCols,resultsDirectory,figureFileName);
        end
        
        if 1
            save([resultsDirectory filesep resultsFileName],'modelResults','modelResultsTable','modelMetaData');
            fprintf('%s\n',['reconXComparisonDriver complete.']);
            fprintf('%s\n',['modelResults saved to ' resultsDirectory filesep resultsFileName]);
        end
end






