% identify ATP contribution(production) rxns
types=fieldnames(Allmodels);
for k=1:length(types)

    %define the samples:
    samples = fieldnames(Allmodels.(types{k}));
    for i=1:length(samples)
        % Perform flux variability analysis and identify the highest contributing reactions
        model=Allmodels.(types{k}).(samples{i});
        
        if FVA==1 
            [minFlux,maxFlux,~, ~] = fastFVA(model);
            %Replace the lower and upper bounds by the min and max Fluxes obtained
            %through FVA:
            model_FVA = model;
            model_FVA.lb = minFlux;
            model_FVA.ub = maxFlux;
        else
            model_FVA=model;
        end
        ResultsAllCellLines.(samples{i}).modelPruned = model_FVA;
        eucNorm = 1e-6;
        %find the contributing rxns for the mets (atp) for each model
        [BMall,ResultsAllCellLines,metRsall,...
            maximum_contributing_rxn,maximum_contributing_flux,...
            ATPyield] = predictFluxSplits_bioenergetics(model_FVA, obj, met2test,samples(i),...
            ResultsAllCellLines, dir,  transportRxns,ATPprod,'EX_glc_D[e]',eucNorm, method);
        if dir == 1
            ATPproduction.(types{k}).(samples{i})=ResultsAllCellLines.(samples{i}).flux_split_atp;
        else
            ATPconsumption.(types{k}).(samples{i})=ResultsAllCellLines.(samples{i}).flux_split_atp;
        end
    end
end