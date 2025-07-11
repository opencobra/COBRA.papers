function [TestSolution,TestSolutionName,TestedRxns,PercTestedRxns] = Test4HumanFctExtv5_eFBA(model,test,optionSinks,param)
%% function [TestSolution,TestSolutionName] = Test4HumanFct(model,test)
% This functions test for the ~460 human functions -  I removed duplicates
%
% INPUT
% model             model structure (Recon1, with desired in silico
%                   condition)
% test              possible statements: Recon1, IECori, IEC, all (default)
%                   (choose IECori if you intend to test the IEC model OR a model that
%                   contains lumen ('u') as compartment otw choose IEC);
%                   all check for Recon1 and IEC
% option            if true = set sink reactions to 0 (default, leave unchanged).
%                   Note that all lb's of exchanges and demands will be set to 0
%
% OUTPUT
% TestSolution      array containing the optimal value for the different
%                   tests
% TestSolutionName  array containing the names  for the different tests
%
% Xi Luo changed 'Test4HumanFctExtv5' function for bioenergeticPD project
%%
if nargin<2
    test = 'all';
end
if nargin<3
    optionSinks = 0; % do not close
end

if nargin<4
    param.debug = 1;
    param.feasTol=1e-7;; % default
end

global saveDiary

if optionSinks
    % close sink reactions
    model.lb(strmatch('sink_',model.rxns))=0;
end

TestSolution = [];
%% Setup

% for organ atlas derived from Harvey only

if saveDiary
    %save each diary to PSCM_toolbox/Files/OrganChecks/
    aPath=which('Test4HumanFctExtv5_eeFBA');
    diary([aPath filesep 'Test4Functions_eeFBA_diary.txt']);
end
TestedRxns =[];
tol = 1e-6;

% skip mets name replacement

% remove all objective function first
model.c(find(model.c)) = 0;
modelOri = model;
k = 1;
RPMI_composition={'EX_ala_L[e]','EX_arg-L[e]','EX_asn_L[e]','EX_asp_L[e]','EX_cys-L[e]','EX_gln-L[e]','EX_glu-L[e]','EX_gly[e]','EX_his-L[e]','EX_ile_L[e]','EX_leu_L[e]','EX_lys-L[e]','EX_met_L[e]','EX_phe_L[e]','EX_4HPRO','EX_pro-L[e]','EX_ser_L[e]','EX_thr_L[e]','EX_trp_L[e]','EX_tyr_L[e]','EX_val_L[e]','EX_ascb_L[e]','EX_btn[e]','EX_chol[e]','EX_pnto_R[e]','EX_fol[e]','EX_ncam[e]','EX_pydxn[e]','EX_ribflv[e]','EX_thm[e]','EX_cbl1[e]','EX_inost[e]','EX_ca2[e]','EX_fe3[e]','EX_k[e]','EX_hco3[e]','EX_na1[e]','EX_pi[e]','EX_glc[e]','EX_hxan[e]','EX_lnlc[e]','EX_lipoate[e]','EX_ptrc[e]','EX_pyr[e]','EX_thymd[e]','EX_etha[e]','EX_gthrd[e]'};

if strcmp(test,'all') 
    % check if model could generate atp under different sources
    %% 1. ATP max aerobic, glc, v0.05
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_glc[e]'))=-1;model.ub(ismember(model.rxns,'EX_glc[e]'))=-1;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    model.c(ismember(model.rxns,'ATPM'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f= sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        if ~isfield(eFBA,'f')
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'ATP max, aerobic, glc';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 2. ATP max, anaerobic glc, v0.05
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_glc[e]'))=-1;model.ub(ismember(model.rxns,'EX_glc[e]'))=-1;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=0;model.ub(ismember(model.rxns,'EX_o2[e]'))=0;
    model.c(ismember(model.rxns,'ATPM'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'ATP max, anaerobic, glc';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 3. ATP max, aerobic, citrate
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_cit[e]'))=-1;model.ub(ismember(model.rxns,'EX_cit[e]'))=-1;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    model.c(ismember(model.rxns,'ATPM'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        if ~isfield(eFBA,'f')
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'ATP max, aerobic, citrate';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 4. ATP max, aerobic, EtOH substrate v0.05
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_etoh[e]'))=-1;model.ub(ismember(model.rxns,'EX_etoh[e]'))=-1;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    model.c(ismember(model.rxns,'ATPM'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        if ~isfield(eFBA,'f')
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'ATP max, aerobic, etoh';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 5. ATP max, aerobic, glutamate v0.05
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_glu-L[e]'))=-1;model.ub(ismember(model.rxns,'EX_glu-L[e]'))=-1;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    model.c(ismember(model.rxns,'ATPM'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        if ~isfield(eFBA,'f')
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'ATP max, aerobic, glu-L';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 6. ATP max, aerobic, glutamine substrate
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_gln-L[e]'))=-1;model.ub(ismember(model.rxns,'EX_gln-L[e]'))=-1;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    model.c(ismember(model.rxns,'ATPM'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        if ~isfield(eFBA,'f')
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'ATP max, aerobic, gln-L';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 7. ATP max, aerobic, glycine substrate v0.05
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_gly[e]'))=-1;model.ub(ismember(model.rxns,'EX_gly[e]'))=-1;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    model.c(ismember(model.rxns,'ATPM'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        if ~isfield(eFBA,'f')
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'ATP max, aerobic, gly';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 8. ATP max, aerobic, lactate substrate v0.05
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_lac-L[e]'))=-1;model.ub(ismember(model.rxns,'EX_lac-L[e]'))=-1;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    model.c(ismember(model.rxns,'ATPM'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        if ~isfield(eFBA,'f')
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'ATP max, aerobic, lac-L';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 9. ATP max, aerobic, proline substrate v0.05
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_pro-L[e]'))=-1;model.ub(ismember(model.rxns,'EX_pro-L[e]'))=-1;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    model.c(ismember(model.rxns,'ATPM'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        if ~isfield(eFBA,'f')
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'ATP max, aerobic, pro-L';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 10. ATP production via default electron transport chain constraints
    model = modelOri;
    model.c(find(model.c)) = 0;
%     model.lb(ismember(model.rxns,'EX_glc[e]'))=-1;model.ub(ismember(model.rxns,'EX_glc[e]'))=0;
    %model.lb(ismember(model.rxns,'CYOOm3'))=1; % there is an alternative
    %reaction
    if ~isempty(find(ismember(model.rxns,'CYOR_u10mi'))) && ~isempty(find(ismember(model.rxns,'NADH2_u10mi')))
        if  model.ub(ismember(model.rxns,'CYOR_u10mi'))>=1 && model.ub(ismember(model.rxns,'NADH2_u10mi'))>=1
            model.c(ismember(model.rxns,'ATPM'))=1;
            if find(model.c)>0
                try
                    eFBA=entropicFluxBalanceAnalysis(model,param);
                    eFBA.f=sum(eFBA.v(find(model.c>0)));
                catch
                    eFBA.f=NaN;
                end
                if ~isfield(eFBA,'f')
                    eFBA.f=NaN;
                end
                TestSolution(k,1) = eFBA.f;
            else
                TestSolution(k,1) = NaN;
            end
        else
            TestSolution(k,1) = NaN;
        end
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'ATP production via electron transport chain';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA   
    %% 11. gthrd reduces h2o2
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_gthrd[e]'))=-1;
    model.c(ismember(model.rxns,'GTHP'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'gthrd reduces h2o2, GTHP [c] ';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    
    model = modelOri;
    model.lb(ismember(model.rxns,'EX_gthrd[e]'))=-1;model.ub(ismember(model.rxns,'gthox[e]'))=1;
    model.c(find(model.c)) = 0;
    model.c(ismember(model.rxns,'GTHPe'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'gthrd reduces h2o2, GTHP [e] ';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_gthrd[e]'))=-1;
    model.c(ismember(model.rxns,'GTHPm'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'gthrd reduces h2o2, GTHP [m] ';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 12. gly -> co2 and nh4 (via glycine cleavage system)
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'gly[c]','co2[c]','nh4[c]'},[-1 -1; 0.1 100; 0.1 100]);
    model.lb(ismember(model.rxns,'EX_nh4[e]'))=0;model.ub(ismember(model.rxns,'EX_nh4[e]'))=1000;
    model.c(ismember(model.rxns,'sink_nh4[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'gly -> co2 + nh4';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    
    %% 13. 12ppd-S -> mthgxl
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'12ppd-S[c]','mthgxl[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_mthgxl[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = '12ppd-S[c] -> mthgxl[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 14. 12ppd-S -> pyr
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'12ppd-S[c]','pyr[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_pyr[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = '12ppd-S[c] -> pyr[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 15. 3pg -> gly
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'3pg[c]','gly[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_gly[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = '3pg[c] -> gly[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 16. 3pg -> ser-L
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'3pg[c]','ser-L[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_ser-L[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = '3pg[c] -> ser-L[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 17. 4abut -> succ[m]
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'4abut[c]','succ[m]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_succ[m]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = '4abut[c] -> succ[m]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 18. 4hpro-LT[m] -> glx[m]
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'4hpro-LT[m]','glx[m]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_glx[m]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = '4hpro-LT[m] -> glx[m]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 19. 5aop -> pheme
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'5aop[c]','pheme[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_pheme[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = '5aop[c] -> pheme[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 20. aact -> mthgxl
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'aact[c]','mthgxl[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_mthgxl[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'aact[c] -> mthgxl[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 21. acac[m] -> acetone[m]
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'acac[m]','acetone[m]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_acetone[m]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'acac[m] -> acetone[m]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 22. acac[m] -> bhb[m]
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'acac[m]','bhb[m]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_bhb[m]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'acac[m] -> bhb[m]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 23. acald -> ac
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'acald[c]','ac[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_ac[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'acald[c] -> ac[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 24. accoa[c] -> pmtcoa[c] -> malcoa[m]
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'accoa[c]','pmtcoa[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_pmtcoa[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'accoa[c] -> pmtcoa[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 25. accoa[c] -> pmtcoa[c] -> malcoa[m]
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'pmtcoa[c]','malcoa[m]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_malcoa[m]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'pmtcoa[c] -> malcoa[m]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 26. acetone -> mthgxl
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'acetone[c]','mthgxl[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_mthgxl[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'acetone[c] -> mthgxl[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 27. acgal -> udpacgal
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'acgal[c]','udpacgal[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_udpacgal[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'acgal[c] -> udpacgal[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 28. acgam -> cmpacna
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'acgam[c]','cmpacna[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_cmpacna[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'acgam[c] -> cmpacna[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 29. acorn -> orn
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'acorn[c]','orn[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_orn[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'acorn[c] -> orn[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 30. adrnl -> 34dhoxpeg
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'adrnl[c]','34dhoxpeg[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_34dhoxpeg[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'adrnl[c] -> 34dhoxpeg[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 31. akg[c] -> glu-L[c] % I adjusted lb since otherwise not feasible
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_akg[e]'))=-1;model.ub(ismember(model.rxns,'EX_akg[e]'))=-1;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    if ~isempty(strmatch('ALATA_L',model.rxns,'exact'))
        model.c(ismember(model.rxns,'ALATA_L'))=1;
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'akg[c] -> glu-L[c] (ALATA_L)';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    
    %% 32. akg[c] -> glu-L[c]
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_akg[e]'))=-1;model.ub(ismember(model.rxns,'EX_akg[e]'))=-1;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    if ~isempty(strmatch('ASPTA',model.rxns,'exact'))
        model.c(ismember(model.rxns,'ASPTA'))=1;
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'akg[c] -> glu-L[c] (ASPTA)';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    
    %% 33. akg[m[ -> oaa[m]
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'akg[m]','oaa[m]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_oaa[m]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'akg[m] -> oaa[m]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 34. akg[m] -> glu-L[m]
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'akg[m]','glu-L[m]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_glu-L[m]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'akg[m] -> glu-L[m]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'akg[m]'},-1 , -1);
    if ~isempty(strmatch('ASPTAm',model.rxns,'exact'))
        model.c(ismember(model.rxns,'ASPTAm'))=1;
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'akg[m] -> glu-L[m] (ASPTAm)';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    
    %% 35. ala-B -> msa
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'ala-B[c]','msa[m]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_msa[m]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'ala-B[c] -> msa[m]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 36. ala-D -> pyr
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'ala-D[c]','pyr[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_pyr[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'ala-D[c] -> pyr[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 37. ala-L -> ala-D
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'ala-L[c]','ala-D[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_ala-D[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'ala-L[c] -> ala-D[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 38. ala-L -> pyr
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'ala-L[c]','pyr[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_pyr[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'ala-L[c] -> pyr[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 39. arachd[c] -> malcoa[m]
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'arachd[c]','malcoa[m]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_malcoa[m]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'arachd[c] -> malcoa[m]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 40. arachd[r] -> txa2[r]
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'arachd[r]','txa2[r]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_txa2[r]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'arachd[r] -> txa2[r]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 41. arg-L -> creat
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'arg-L[c]','creat[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_creat[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'arg-L[c] -> creat[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 42. arg-L -> glu-L [m]
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'arg-L[c]','glu-L[m]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_glu-L[m]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'arg-L -> glu-L [m]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 43. arg-L -> no
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'arg-L[c]','no[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_no[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'arg-L -> no';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 44. arg-L -> pcreat
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'arg-L[c]','pcreat[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_pcreat[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'arg-L[c] -> pcreat[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 45. ascb -> eryth
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'ascb-L[c]','eryth[c]'},[-1 -1; 0 100]);
    model.lb(find(ismember(model.rxns,'DM_ascb_L[c]')))=-1;
    model.ub(find(ismember(model.rxns,'DM_ascb_L[c]')))=-1;
    model.c(ismember(model.rxns,'sink_eryth[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'ascb-L[c] -> eryth[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 46. ascb -> lyxnt
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'ascb-L[c]','lyxnt[c]'},[-1 -1; 0 100]);
    model.lb(find(ismember(model.rxns,'DM_ascb_L[c]')))=-1;
    model.ub(find(ismember(model.rxns,'DM_ascb_L[c]')))=-1;
    model.c(ismember(model.rxns,'sink_lyxnt[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'ascb-L[c] -> lyxnt[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 47. ascb -> thrnt
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'ascb-L[c]','thrnt[c]'},[-1 -1; 0 100]);
    model.lb(find(ismember(model.rxns,'DM_ascb_L[c]')))=-1;
    model.ub(find(ismember(model.rxns,'DM_ascb_L[c]')))=-1;
    model.c(ismember(model.rxns,'sink_thrnt[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'ascb-L[c] -> thrnt[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 48. ascb -> xylnt
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'ascb-L[c]','xylnt[c]'},[-1 -1; 0 100]);
    model.lb(find(ismember(model.rxns,'DM_ascb_L[c]')))=-1;
    model.ub(find(ismember(model.rxns,'DM_ascb_L[c]')))=-1;
    model.c(ismember(model.rxns,'sink_xylnt[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'ascb-L[c] -> xylnt[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 49. asn-L -> oaa
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'asn-L[c]','oaa[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_oaa[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'asn-L[c] -> oaa[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 50. asp-L + hco3 -> arg-L
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'asp-L[c]','hco3[c]','arg-L[c]'},[-1 -1;-1 -1;0 100]);
    model.c(ismember(model.rxns,'sink_arg-L[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'asp-L[c] + hco3[c] -> arg-L[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 51. asp-L -> ala-B
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'asp-L[c]','ala-B[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_ala-B[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'asp-L[c] -> ala-B[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 52. asp-L -> asn-L
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'asp-L[c]','asn-L[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_asn-L[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'asp-L[c] -> asn-L[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 53. asp-L -> fum (via argsuc)
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'asp-L[c]','argsuc[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_argsuc[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'asp-L[c] -> argsuc[c], asp-L -> fum (via argsuc), 1';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'argsuc[c]','fum[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_fum[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'argsuc[c] -> fum[c], asp-L -> fum (via argsuc), 2';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 54. asp-L -> fum (via dcamp)
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'asp-L[c]','dcamp[c]'},[-1 -1; 0 100]);
    model.lb(find(ismember(model.rxns,'sink_asp_L[c]')))=-1;
    model.ub(find(ismember(model.rxns,'sink_asp_L[c]')))=-1;
    model.c(ismember(model.rxns,'sink_dcamp[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'asp-L[c] -> dcamp[c], asp-L -> fum (via dcamp), 1';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(find(ismember(model.rxns,'sink_asp_L[c]')))=-1;
    model.ub(find(ismember(model.rxns,'sink_asp_L[c]')))=-1;
    [model] = addSinkReactions(model,{'dcamp[c]','fum[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_fum[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'dcamp[c] -> fum[c], asp-L -> fum (via dcamp), 2';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(find(ismember(model.rxns,'sink_asp_L[c]')))=-1;
    model.ub(find(ismember(model.rxns,'sink_asp_L[c]')))=-1;
    [model] = addSinkReactions(model,{'dcamp[c]','fum[c]'},[-1 -1; 0 100]);
    if ~isempty(strmatch('ADSS',model.rxns,'exact'))
        model.c(ismember(model.rxns,'ADSS'))=1;
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'dcamp[c] -> fum[c], asp-L -> fum (via dcamp), 3';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    
    %% 55. asp-L -> oaa
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'asp-L[c]','oaa[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_oaa[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'asp-L[c] -> oaa[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 56. carn -> ala-B
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'carn[c]','ala-B[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_ala-B[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'carn -> ala-B';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 57. chol[c] + dag_hs[c] -> pe_hs[c]
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'chol[c]','dag_hs[c]','pe_hs[c]'},[-1 -1;-1 -1;0 100]);
    model.c(ismember(model.rxns,'sink_pe_hs[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'chol[c] + dag_hs[c] -> pe_hs[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 58. choline -> betaine -> glycine
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'chol[m]','glyb[m]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_glyb[m]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'choline -> betaine (glyb) -> glycine, 1 [m]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'glyb[m]','gly[m]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_gly[m]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'choline -> betaine (glyb) -> glycine, 2 [m]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 59. coke[r] -> pecgoncoa[r]
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_glc[e]'))=-1;model.ub(ismember(model.rxns,'EX_glc[e]'))=0;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    [model] = addSinkReactions(model,{'coke[r]','pecgoncoa[r]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_pecgoncoa[r]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'coke[r] -> pecgoncoa[r]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 60. core2[g] -> ksii_core2[g]
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'core2[g]','ksii_core2[g]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_ksii_core2[g]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'core2[g] -> ksii_core2[g]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 61. core4[g] -> ksii_core4[g]
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'core4[g]','ksii_core4[g]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_ksii_core4[g]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'core4[g] -> ksii_core4[g]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 62. cspg_a[ly] -> 2 gal[ly] + glcur[ly] + xyl-D[ly] %I adjusted lb since otw infeasible
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'cspg_a[l]','gal[l]','glcur[l]','xyl-D[l]'},[-1 -1; 0.1 100; 0.1 100; 0.1 100]);
    model.c(ismember(model.rxns,'sink_xyl-D[l]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'cspg_a[ly] -> gal[ly] + glcur[ly] + xyl-D[ly]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 63. cspg_b[ly] -> 2gal[ly] + glcur[ly] + xyl-D[ly]
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'cspg_b[l]','gal[l]','glcur[l]','xyl-D[l]'},[-1 -1; 0.1 100; 0.1 100; 0.1 100]);
    model.c(ismember(model.rxns,'sink_xyl-D[l]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'cspg_b[ly] -> gal[ly] + glcur[ly] + xyl-D[ly]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 64. cspg_c[ly] -> 2 gal[ly] + glcur[ly] + xyl-D[ly]
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'cspg_c[l]','gal[l]','glcur[l]','xyl-D[l]'},[-1 -1; 0.1 100; 0.1 100; 0.1 100]);
    model.c(ismember(model.rxns,'sink_xyl-D[l]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'cspg_c[ly] -> gal[ly] + glcur[ly] + xyl-D[ly]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 65. cspg_d[ly] -> 2 gal[ly] + glcur[ly] + xyl-D[ly]
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'cspg_d[l]','gal[l]','glcur[l]','xyl-D[l]'},[-1 -1; 0.1 100; 0.1 100; 0.1 100]);
    model.c(ismember(model.rxns,'sink_xyl-D[l]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'cspg_d[ly] -> gal[ly] + glcur[ly] + xyl-D[ly]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 66. cspg_e[ly] -> 2 gal[ly] + glcur[ly] + xyl-D[ly]
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'cspg_e[l]','gal[l]','glcur[l]','xyl-D[l]'},[-1 -1; 0.1 100; 0.1 100; 0.1 100]);
    model.c(ismember(model.rxns,'sink_xyl-D[l]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'cspg_e[ly] -> gal[ly] + glcur[ly] + xyl-D[ly]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 67. cys-L + glu-L + gly -> ghtrd
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'cys-L[c]', 'glu-L[c]','gly[c]','gthrd[c]'},[-1 -1;-1 -1;-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_gthrd[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'cys-L + glu-L + gly -> ghtrd';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 68. cys-L -> 3sala -> so4 %I adjusted lb since otw infeasible
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'cys-L[c]','3sala[c]'},[-1 -1; 0 100]);
    model.lb(ismember(model.rxns,'EX_so4[e]'))=0;model.ub(ismember(model.rxns,'EX_so4[e]'))=1000;
    model.c(ismember(model.rxns,'sink_3sala[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'cys-L -> 3sala -> so4, 1';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'3sala[c]','so4[c]'},[-1 -1; 0 100]);
    model.lb(ismember(model.rxns,'EX_so4[e]'))=0;model.ub(ismember(model.rxns,'EX_so4[e]'))=1000;
    model.c(ismember(model.rxns,'sink_so4[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'cys-L -> 3sala -> so4, 2';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 69. cys-L -> hyptaur
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'cys-L[c]','hyptaur[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_hyptaur[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'cys-L[c] -> hyptaur[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 70. cystine -> cys-L
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'Lcystin[c]','cys-L[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_cys-L[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'cystine (Lcystin) -> cys-L';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 71. dhap -> mthgxl
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'dhap[c]','mthgxl[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_mthgxl[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'dhap[c] -> mthgxl[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 72. dmpp -> ggdp
    model = modelOri;
    model.c(find(model.c)) = 0;
    
    for i = 1 : length(RPMI_composition)
        model = changeRxnBounds(model,RPMI_composition{i},-1,'l');
    end
    
    [model] = addSinkReactions(model,{'dmpp[c]','ggdp[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_ggdp[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'dmpp[c] -> ggdp[c] (with RPMI medium)';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 73. dna(n) -> dna5mtc(n)
    model = modelOri;
    model.c(find(model.c)) = 0;
    for i = 1 : length(RPMI_composition)
        model = changeRxnBounds(model,RPMI_composition{i},-1,'l');
    end
    [model] = addSinkReactions(model,{'dna(n)','dna5mtc(n)'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_dna5mtc(n)'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'dna(n) -> dna5mtc(n) (with RPMI medium)';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 74. dolichol_L -> dolmanp_L[r]
    model = modelOri;
    model.c(find(model.c)) = 0;
    for i = 1 : length(RPMI_composition)
        model = changeRxnBounds(model,RPMI_composition{i},-1,'l');
    end
    [model] = addSinkReactions(model,{'dolichol_L[c]','dolmanp_L[r]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_dolmanp_L[r]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'dolichol_L[c] -> dolmanp_L[r] (with RPMI medium)';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 75. dolichol_L -> g3m8mpdol_L[r]
    model = modelOri;
    model.c(find(model.c)) = 0;
    for i = 1 : length(RPMI_composition)
        model = changeRxnBounds(model,RPMI_composition{i},-1,'l');
    end
    [model] = addSinkReactions(model,{'dolichol_L[c]','g3m8mpdol_L[r]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_g3m8mpdol_L[r]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'dolichol_L[c] -> g3m8mpdol_L[r] (with RPMI medium)';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 76. dolichol_U -> dolmanp_U[r]
    model = modelOri;
    for i = 1 : length(RPMI_composition)
        model = changeRxnBounds(model,RPMI_composition{i},-1,'l');
    end
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'dolichol_U[c]','dolmanp_U[r]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_dolmanp_U[r]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'dolichol_U[c] -> dolmanp_U[r] (with RPMI medium)';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 77. dolichol_U -> g3m8mpdol_U[r]
    model = modelOri;
    for i = 1 : length(RPMI_composition)
        model = changeRxnBounds(model,RPMI_composition{i},-1,'l');
    end
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'dolichol_U[c]','g3m8mpdol_U[r]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_g3m8mpdol_U[r]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'dolichol_U[c] -> g3m8mpdol_U[r] (with RPMI medium)';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 78. dopa -> homoval (1)
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_glc[e]'))=-1;model.ub(ismember(model.rxns,'EX_glc[e]'))=0;
    [model] = addSinkReactions(model,{'dopa[c]','homoval[c]'},[-1 -1; 0 100]);
    model.lb(find(ismember(model.rxns,'DM_dopa[c]')))=-1;
    model.ub(find(ismember(model.rxns,'DM_dopa[c]')))=-1;
    model.c(ismember(model.rxns,'sink_homoval[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'dopa[c] -> homoval[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 79. etoh -> acald
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'etoh[c]','acald[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_acald[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'etoh[c] -> acald[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 80. f6p + g3p -> r5p
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'f6p[c]','g3p[c]','r5p[c]'},[-1 -1; -1 -1;0 100]);
    model.c(ismember(model.rxns,'sink_r5p[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'f6p[c] + g3p[c] -> r5p[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 81. frdp -> dolichol_L
    model = modelOri;
    model.c(find(model.c)) = 0;
    for i = 1 : length(RPMI_composition)
        model = changeRxnBounds(model,RPMI_composition{i},-1,'l');
    end
    [model] = addSinkReactions(model,{'frdp[c]','dolichol_L[r]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_dolichol_L[r]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'frdp[c] -> dolichol_L[r] (with RPMI medium)';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 82. frdp -> dolichol_U
    model = modelOri;
    model.c(find(model.c)) = 0;
    for i = 1 : length(RPMI_composition)
        model = changeRxnBounds(model,RPMI_composition{i},-1,'l');
    end
    [model] = addSinkReactions(model,{'frdp[c]','dolichol_U[r]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_dolichol_U[r]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'frdp[c] -> dolichol_U[r] (with RPMI medium)';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 83. from ade[c] to amp[c]
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_glc[e]'))=-1;model.ub(ismember(model.rxns,'EX_glc[e]'))=0;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    [model] = addSinkReactions(model,{'ade[c]','amp[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_amp[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'ade[c] -> amp[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 84. from adn[c] to urate[x]
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'adn[c]','urate[x]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_urate[x]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'adn[c] -> urate[x]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 85. from ADP[c] to dATP(n)
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model,rxnsInModel] = addSinkReactions(model,{'adp[c]','datp(n)'},[-1 -1; 0 100]);
    if (rxnsInModel(2) >-1) % reaction exits already in model
        model=changeObjective(model,model.rxns(rxnsInModel(2),1));
    else
        model=changeObjective(model,'sink_datp(n)',1);
    end
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'adp[c] -> datp(n)';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 86. from CDP[c] to dCTP(n)
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model,rxnsInModel] = addSinkReactions(model,{'cdp[c]','dctp(n)'},[-1 -1; 0 100]);
    if (rxnsInModel(2) >-1) % reaction exits already in model
        model=changeObjective(model,model.rxns(rxnsInModel(2),1));
    else
        model=changeObjective(model,'sink_dctp(n)',1);
    end
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'cdp[c] -> dctp(n)';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 87. from cmp to cytd
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'cmp[c]','cytd[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_cytd[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'cmp[c] -> cytd[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 88. from cytd to ala-B
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'cytd[c]','ala-B[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_ala-B[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'cytd[c] -> ala-B[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 89. from dcmp to ala-B
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'dcmp[c]','ala-B[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_ala-B[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'dcmp[c] -> ala-B[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 90. from GDP[c] to dGTP(n)
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'gdp[c]','dgtp(n)'},[-1 -1; 0 100]);
    if (rxnsInModel(2) >-1) % reaction exits already in model
        model=changeObjective(model,model.rxns(rxnsInModel(2),1));
    else
        model=changeObjective(model,'sink_dgtp(n)',1);
    end
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'gdp[c] -> dgtp(n)';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 91. from gln-L + HCO3 to UMP[c]
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'gln-L[c]','hco3[c]','ump[c]'},[-1 -1;-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_ump[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'gln-L + HCO3 -> UMP[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 92. from gsn[c] to urate[x]
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'gsn[c]','urate[x]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_urate[x]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'gsn[c] -> urate[x]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 93. from gua[c] to gmp[c]
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_glc[e]'))=-1;model.ub(ismember(model.rxns,'EX_glc[e]'))=0;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    [model] = addSinkReactions(model,{'gua[c]','gmp[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_gmp[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'gua[c] -> gmp[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 94. from hxan[c] to imp[c]
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_glc[e]'))=-1;model.ub(ismember(model.rxns,'EX_glc[e]'))=0;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    [model] = addSinkReactions(model,{'hxan[c]','imp[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_imp[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'hxan[c] -> imp[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 95. from imp to ATP
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'imp[c]','atp[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_atp[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'imp[c] -> atp[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 96. from imp to gtp
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'imp[c]','gtp[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_gtp[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'imp[c] -> gtp[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 97. from imp[c] to urate[x]
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'imp[c]','urate[x]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_urate[x]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'imp[c] -> urate[x]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 98. from prpp to imp
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'prpp[c]','imp[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_imp[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'prpp[c] -> imp[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 99. from pydx[c] to pydx5p[c]
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_glc[e]'))=-1;model.ub(ismember(model.rxns,'EX_glc[e]'))=0;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    [model] = addSinkReactions(model,{'pydx[c]','pydx5p[c]'},[-1 -1; 0 100]);
    model.lb(find(ismember(model.rxns,'sink_pydx[c]')))=-1;
    model.ub(find(ismember(model.rxns,'sink_pydx[c]')))=-1;
    model.c(ismember(model.rxns,'sink_pydx5p[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'pydx[c] -> pydx5p[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 100. from thm[c] to thmpp[c]
    model = modelOri;
    model.c(find(model.c)) = 0;
    for i = 1 : length(RPMI_composition)
        model = changeRxnBounds(model,RPMI_composition{i},-1,'l');
    end
    [model] = addSinkReactions(model,{'thm[c]','thmpp[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_thmpp[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'thm[c] -> thmpp[c] (with RPMI medium)';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 101. from thm[e] to thmpp[m] %does not work; changing lb has no effect
    model = modelOri;
    model.c(find(model.c)) = 0;
    for i = 1 : length(RPMI_composition)
        model = changeRxnBounds(model,RPMI_composition{i},-1,'l');
    end
    [model] = addSinkReactions(model,{'thm[e]','thmpp[m]'},[-1 -1; 0 100]);
    model.lb(find(ismember(model.rxns,'EX_thm[e]')))=-1;
    model.ub(find(ismember(model.rxns,'EX_thm[e]')))=-1;
    model.c(ismember(model.rxns,'sink_thmpp[m]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'thm[e] -> thmpp[m] (with RPMI medium)';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 102. from thmmp[e] to thmpp[c]
    model = modelOri;
    model.c(find(model.c)) = 0;
    for i = 1 : length(RPMI_composition)
        model = changeRxnBounds(model,RPMI_composition{i},-1,'l');
    end
    [model] = addSinkReactions(model,{'thmmp[e]','thmpp[c]'},[-1 -1; 0 100]);
    model.lb(find(ismember(model.rxns,'EX_thmmp[e]')))=-1;
    model.ub(find(ismember(model.rxns,'EX_thmmp[e]')))=-1;
    model.c(ismember(model.rxns,'sink_thmpp[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'thmmp[e] -> thmpp[c] (with RPMI medium)';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 103. from thmmp[e] to thmpp[m]
    model = modelOri;
    model.c(find(model.c)) = 0;
    for i = 1 : length(RPMI_composition)
        model = changeRxnBounds(model,RPMI_composition{i},-1,'l');
    end
    model.lb(find(ismember(model.rxns,'EX_thmmp[e]')))=-1;
    model.ub(find(ismember(model.rxns,'EX_thmmp[e]')))=-1;
    [model] = addSinkReactions(model,{'thmpp[m]'},[ 0 100]);
    model.c(ismember(model.rxns,'sink_thmpp[m]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'thmmp[e] -> thmpp[m] (with RPMI medium)';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 104. from tyr-L[m] to q10[m]
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'tyr-L[m]','q10[m]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_q10[m]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'tyr-L[m] -> q10[m]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 105. from UDP[c] to dTTP(n)
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'udp[c]','dttp(n)'},[-1 -1; 0 100]);
    if (rxnsInModel(2) >-1) % reaction exits already in model
        model=changeObjective(model,model.rxns(rxnsInModel(2),1));
    else
        model=changeObjective(model,'sink_dttp(n)',1);
    end
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'udp[c] -> dttp(n)';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 106. from ump to ala-B
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'ump[c]','ala-B[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_ala-B[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'ump[c] -> ala-B[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 107. fru -> dhap
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'fru[c]','dhap[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_dhap[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'fru[c] -> dhap[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 108. fru -> g3p
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'fru[c]','g3p[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_g3p[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'fru[c] -> g3p[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 109. fuc -> gdpfuc
    model = modelOri;
    model.c(find(model.c)) = 0;
    
    for i = 1 : length(RPMI_composition)
        model = changeRxnBounds(model,RPMI_composition{i},-1,'l');
    end
    [model] = addSinkReactions(model,{'fuc-L[c]','gdpfuc[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_gdpfuc[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'fuc-L[c] -> gdpfuc[c] (with RPMI medium)';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 110. fum[m] -> oaa[m]
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'fum[m]','oaa[m]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_oaa[m]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'fum[m] -> oaa[m]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 111. g1p -> dtdprmn
    model = modelOri;
    model.c(find(model.c)) = 0;
    for i = 1 : length(RPMI_composition)
        model = changeRxnBounds(model,RPMI_composition{i},-1,'l');
    end
    [model] = addSinkReactions(model,{'g1p[c]','dtdprmn[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_dtdprmn[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'g1p[c] -> dtdprmn[c] (with RPMI medium)';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 112. g3p -> mthgxl
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'g3p[c]','mthgxl[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_mthgxl[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'g3p[c] -> mthgxl[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 113. g6p -> r5p
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'g6p[c]','r5p[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_r5p[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'g6p[c] -> r5p[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 114. g6p -> ru5p
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'g6p[c]','ru5p-D[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_ru5p-D[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'g6p[c] -> ru5p-D[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 115. gal -> glc
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'gal[c]','glc-D[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_glc-D[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'gal[c] -> glc-D[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 116. gal -> udpgal
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'gal[c]','udpgal[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_udpgal[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'gal[c] -> udpgal[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 117. galgluside[g] -> galgalgalthcrm_hs[g]
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'galgluside_hs[g]','galgalgalthcrm_hs[g]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_galgalgalthcrm_hs[g]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'galgluside_hs[g] -> galgalgalthcrm_hs[g]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 118. galgluside_hs[g] -> acgagbside_hs[g]
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'galgluside_hs[g]','acgagbside_hs[g]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_acgagbside_hs[g]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'galgluside_hs[g] -> acgagbside_hs[g]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 119. galgluside_hs[g] -> acnacngalgbside_hs[g]
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'galgluside_hs[g]','acnacngalgbside_hs[g]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_acnacngalgbside_hs[g]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'galgluside_hs[g] -> acnacngalgbside_hs[g]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 120. galgluside_hs[g] -> gd1b2_hs[g]
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'galgluside_hs[g]','gd1b2_hs[g]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_gd1b2_hs[g]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'galgluside_hs[g] -> gd1b2_hs[g]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 121. galgluside_hs[g] -> gd1c_hs[g]
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'galgluside_hs[g]','gd1c_hs[g]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_gd1c_hs[g]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'galgluside_hs[g] -> gd1c_hs[g]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 122. galgluside_hs[g] -> gp1c_hs[g]
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'galgluside_hs[g]','gp1c_hs[g]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_gp1c_hs[g]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'galgluside_hs[g] -> gp1c_hs[g]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 123. galgluside_hs[g] -> gq1balpha_hs[g]
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'galgluside_hs[g]','gq1balpha_hs[g]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_gq1balpha_hs[g]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'galgluside_hs[g] -> gq1balpha_hs[g]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 124. gam6p -> uacgam
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'gam6p[c]','uacgam[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_uacgam[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'gam6p[c] -> uacgam[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 125. gdpmann -> gdpfuc
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_glc[e]'))=-1;model.ub(ismember(model.rxns,'EX_glc[e]'))=0;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    [model] = addSinkReactions(model,{'gdpmann[c]','gdpfuc[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_gdpfuc[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'gdpmann[c] -> gdpfuc[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 126. glc -> inost
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'glc-D[c]','inost[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_inost[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'glc-D[c] -> inost[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 127. glc -> lac + atp + h2o % I assumed lac-L
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'glc-D[c]','lac-L[c]','atp[c]','h2o[c]'},[-1 -1; 0.1 100; 0.1 100; 0.1 100]);
    model.c(ismember(model.rxns,'sink_lac-L[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'glc-D[c] -> lac-L[c] + atp[c] + h2o[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 128. glc -> lac-D
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'glc-D[c]','lac-D[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_lac-D[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'glc-D[c] -> lac-D[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 129. glc -> lcts[g] (2)
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'glc-D[c]','lcts[g]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_lcts[g]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'glc-D[c] -> lcts[g]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 130. glc -> pyr
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'glc-D[c]','pyr[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_pyr[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'glc-D[c] -> pyr[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 131. gln -> nh4
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'gln-L[c]','nh4[c]'},[-1 -1; 0 100]);
    model.lb(ismember(model.rxns,'EX_nh4[e]'))=0;model.ub(ismember(model.rxns,'EX_nh4[e]'))=1000;
    model.c(ismember(model.rxns,'sink_nh4[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'gln-L[c] -> nh4[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 132. gln-L[m] -> glu-L[m]
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'gln-L[m]','glu-L[m]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_glu-L[m]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'gln-L[m] -> glu-L[m]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 133. glu5sa -> pro-L
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'glu5sa[c]','pro-L[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_pro-L[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'glu5sa[c] -> pro-L[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 134. glu-L -> 4abut
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'glu-L[c]','4abut[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_4abut[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'glu-L[c] -> 4abut[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 135. glu-L -> gln-L[c]
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'glu-L[c]','gln-L[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_gln-L[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'glu-L[c] -> gln-L[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 136. glu-L -> pro-L
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'glu-L[c]','pro-L[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_pro-L[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'glu-L -> pro-L';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 137. glu-L[m] -> akg[m]
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'glu-L[m]','akg[m]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_akg[m]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'glu-L[m] -> akg[m]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 138. gluside_hs[g] -> galgluside_hs[g]
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'gluside_hs[g]','galgluside_hs[g]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_galgluside_hs[g]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'gluside_hs[g] -> galgluside_hs[g]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 139. glx[m] -> glyclt[m]
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'glx[m]','glyclt[m]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_glyclt[m]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'glx[m] -> glyclt[m]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 140. gly -> ser-L -> pyr (via SERD_L) % SERD_L does not exist in human
    % (L-serine deaminase)
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'gly[c]','ser-L[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_ser-L[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'gly[c] -> ser-L[c] -> pyr[c], 1';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'ser-L[c]','pyr[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_pyr[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'gly[c] -> ser-L[c] -> pyr[c], 2';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 141. glyc -> glc
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'glyc[c]','glc-D[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_glc-D[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'glyc[c] -> glc-D[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 142. glyc[c] + Rtotal[c] + Rtotal2[c] -> dag_hs[c]
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'glyc[c]','Rtotal[c]','Rtotal2[c]','dag_hs[c]'},[-1 -1;-1 -1;-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_dag_hs[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'glyc[c] + Rtotal[c] + Rtotal2[c] -> dag_hs[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 143. glyc[c] + Rtotal[c] -> tag_hs[c]
    model = modelOri;
    model.c(find(model.c)) = 0;
    for i = 1 : length(RPMI_composition)
        model = changeRxnBounds(model,RPMI_composition{i},-1,'l');
    end
    [model] = addSinkReactions(model,{'glyc[c]','Rtotal[c]','tag_hs[c]'},[-1 -1;-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_tag_hs[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'glyc[c] + Rtotal[c] -> tag_hs[c] (with RPMI medium)';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 144. glyclt -> gly
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'glyclt[c]','gly[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_gly[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'glyclt[c] -> gly[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 145. glygn2 -> glc % changing lb has no effect
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'glygn2[c]','glc-D[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_glc-D[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'glygn2[c] -> glc-D[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 146. glygn2[e] -> glc[e]
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'glygn2[e]','glc-D[e]'},[-1 -1; 0 100]);
    %model.c(ismember(model.rxns,'sink_glc-D[e]'))=1;
    if ~isempty(strmatch('AMY2e',model.rxns,'exact'))
        model = changeObjective(model,'AMY2e',1);
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'glygn2[e] -> glc-D[e] - via AMY2e';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 147. glyx -> oxa % I assumed glx
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'glx[c]','oxa[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_oxa[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'glx[c] -> oxa[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 148. ha[l] -> acgam[l] + glcur[l]
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'ha[l]','acgam[l]','glcur[l]'},[-1 -1; 0.1 100; 0.1 100]);
    model.c(ismember(model.rxns,'sink_acgam[l]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'ha[l] -> acgam[l] + glcur[l]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 149. His -> glu-L
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'his-L[c]','glu-L[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_glu-L[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'his-L[c] -> glu-L[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 150. his-L -> hista
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'his-L[c]','hista[c]'},[-1 -1; 0 100]);
    model.lb(find(ismember(model.rxns,'sink_his_L[c]')))=-1;
    model.ub(find(ismember(model.rxns,'sink_his_L[c]')))=-1;
    model.c(ismember(model.rxns,'sink_hista[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'his-L[c] -> hista[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 151. hista -> 3mlda
    model = modelOri;
    model.c(find(model.c)) = 0;
    for i = 1 : length(RPMI_composition)
        model = changeRxnBounds(model,RPMI_composition{i},-1,'l');
    end
    [model] = addSinkReactions(model,{'hista[c]','3mlda[c]'},[-1 -1; 0 100]);
    model.lb(find(ismember(model.rxns,'DM_hista[c]')))=-1;
    model.ub(find(ismember(model.rxns,'DM_hista[c]')))=-1;
    model.c(ismember(model.rxns,'sink_3mlda[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'hista[c] -> 3mlda[c] (with RPMI medium)';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 152. hista -> im4ac
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'hista[c]','im4act[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_im4act[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'hista[c] -> im4ac[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 153. hmgcoa[x] -> chsterol[r]
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'hmgcoa[x]','chsterol[r]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_chsterol[r]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'hmgcoa[x] -> chsterol[r]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 154. hmgcoa[x] -> frdp[x]
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'hmgcoa[x]','frdp[x]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_frdp[x]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'hmgcoa[x] -> frdp[x]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 155. hmgcoa[x] -> xoldiolone[r]
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'hmgcoa[x]','xoldiolone[r]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_xoldiolone[r]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'hmgcoa[x] -> xoldiolone[r]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 156. hmgcoa[x] -> xoltriol[c]
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'hmgcoa[x]','xoltriol[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_xoltriol[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'hmgcoa[x] -> xoltriol[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 157. hpyr -> 2pg
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'hpyr[c]','2pg[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_2pg[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'hpyr[c] -> 2pg[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 158. hpyr -> glyclt
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'hpyr[c]','glyclt[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_glyclt[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'hpyr[c] -> glyclt[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 159. hpyr -> glyc-S
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'hpyr[c]','glyc-S[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_glyc-S[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'hpyr[c] -> glyc-S[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 160. hspg[l] -> 2 gal[l] + glcur[l] + xyl-D[l]
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'hspg[l]','gal[l]','glcur[l]','xyl-D[l]'},[-1 -1; 0.1 100; 0.1 100; 0.1 100]);
    model.c(ismember(model.rxns,'sink_xyl-D[l]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'hspg[l] -> gal[l] + glcur[l] + xyl-D[l]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 161. hyptaur[c] -> taur[x]
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'hyptaur[c]','taur[x]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_taur[x]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'hyptaur[c] -> taur[x]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 162. ile-L -> accoa
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'ile-L[c]','accoa[c]'},[-1 -1; 0 100]);
    model.lb(find(ismember(model.rxns,'sink_ile_L[c]')))=-1;
    model.ub(find(ismember(model.rxns,'sink_ile_L[c]')))=-1;
    model.lb(find(ismember(model.rxns,'sink_coa[c]')))=-1;
    model.ub(find(ismember(model.rxns,'sink_coa[c]')))=1;
    model.c(ismember(model.rxns,'sink_accoa[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'ile-L[c] -> accoa[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 163. inost -> pail_hs
    model = modelOri;
    model.c(find(model.c)) = 0;
    for i = 1 : length(RPMI_composition)
        model = changeRxnBounds(model,RPMI_composition{i},-1,'l');
    end
    [model] = addSinkReactions(model,{'inost[c]','pail_hs[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_pail_hs[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'inost[c] -> pail_hs[c] (with RPMI medium)';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 164. inost -> pail45p_hs
    model = modelOri;
    model.c(find(model.c)) = 0;
    for i = 1 : length(RPMI_composition)
        model = changeRxnBounds(model,RPMI_composition{i},-1,'l');
    end
    [model] = addSinkReactions(model,{'inost[c]','pail45p_hs[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_pail45p_hs[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'inost[c] -> pail45p_hs[c] (with RPMI medium)';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 165. inost -> pail4p_hs
    model = modelOri;
    for i = 1 : length(RPMI_composition)
        model = changeRxnBounds(model,RPMI_composition{i},-1,'l');
    end
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'inost[c]','pail4p_hs[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_pail4p_hs[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'inost[c] -> pail4p_hs[c] (with RPMI medium)';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 166. inost -> xu5p-D
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'inost[c]','xu5p-D[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_xu5p-D[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'inost[c] -> xu5p-D[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 167. ipdp[x] -> sql[r]
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_glc[e]'))=-1;model.ub(ismember(model.rxns,'EX_glc[e]'))=0;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    [model] = addSinkReactions(model,{'ipdp[x]','sql[r]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_sql[r]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'ipdp[x] -> sql[r]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 168. itacon[m] -> pyr[m] %changing lb has no effect
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'itacon[m]','pyr[m]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_pyr[m]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'itacon[m] -> pyr[m]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 169. ksi[l] -> man[l] + acgam[l]
    model = modelOri;
    model.c(find(model.c)) = 0;
    for i = 1 : length(RPMI_composition)
        model = changeRxnBounds(model,RPMI_composition{i},-1,'l');
    end
    [model] = addSinkReactions(model,{'ksi[l]','man[l]','acgam[l]'},[-1 -1; 0.1 100; 0.1 100]);
    model.c(ismember(model.rxns,'sink_acgam[l]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'ksi[l] -> man[l] + acgam[l] (with RPMI medium)';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 170. ksii_core2[l] -> Ser/Thr[l]
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model,rxnsInModel] = addSinkReactions(model,{'ksii_core2[l]','Ser/Thr[l]'},[-1 -1; 0 100]);
    if (rxnsInModel(2) >-1) % reaction exits already in model
        model=changeObjective(model,model.rxns(rxnsInModel(2),1));
    else
        model=changeObjective(model,'sink_Ser/Thr[l]',1);
    end
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'ksii_core2[l] -> Ser/Thr[l]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 171. ksii_core4[l] -> Ser/Thr[l]
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model,rxnsInModel] = addSinkReactions(model,{'ksii_core4[l]','Ser/Thr[l]'},[-1 -1; 0 100]);
    if (rxnsInModel(2) >-1) % reaction exits already in model
        model=changeObjective(model,model.rxns(rxnsInModel(2),1));
    else
        model=changeObjective(model,'sink_Ser/Thr[l]',1);
    end
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'ksii_core4[l] -> Ser/Thr[l]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 172. l2fn2m2masn[g] -> ksi[g]
    model = modelOri;
    model.c(find(model.c)) = 0;
    for i = 1 : length(RPMI_composition)
        model = changeRxnBounds(model,RPMI_composition{i},-1,'l');
    end
    [model] = addSinkReactions(model,{'l2fn2m2masn[g]','ksi[g]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_ksi[g]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'l2fn2m2masn[g] -> ksi[g] (with RPMI medium)';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 173. lac -> glc % i assumed lac-L
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'lac-L[c]','glc-D[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_glc-D[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'lac-L[c] -> glc-D[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 174. Lcyst[c] -> taur[x]
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'Lcyst[c]','taur[x]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_taur[x]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Lcyst[c] -> taur[x]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 175. leu-L -> accoa
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'leu-L[c]','accoa[c]'},[-1 -1; 0 100]);
    model.lb(find(ismember(model.rxns,'sink_leu_L[c]')))=-1;
    model.ub(find(ismember(model.rxns,'sink_leu_L[c]')))=-1;
    model.lb(find(ismember(model.rxns,'sink_coa[c]')))=-1;
    model.ub(find(ismember(model.rxns,'sink_coa[c]')))=1;
    model.c(ismember(model.rxns,'sink_accoa[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'leu-L[c] -> accoa[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 176. lys-L[c] -> accoa[m] (via saccrp-L pathway)
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'lys-L[c]','accoa[m]'},[-1 -1; 0 100]);
    model.lb(find(ismember(model.rxns,'sink_lys_L[c]')))=-1;
    model.ub(find(ismember(model.rxns,'sink_lys_L[c]')))=-1;
    model.lb(find(ismember(model.rxns,'sink_coa[c]')))=-1;
    model.ub(find(ismember(model.rxns,'sink_coa[c]')))=1;
    model.c(ismember(model.rxns,'sink_accoa[m]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'lys-L[c] -> accoa[m] (via saccrp-L pathway)';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 177. lys-L[x] -> aacoa[m] (via Lpipecol pathway)
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(find(ismember(model.rxns,'EX_lys-L[e]')))=-1;
    model.ub(find(ismember(model.rxns,'EX_lys-L[e]')))=-1;
    [model] = addSinkReactions(model,{'lys-L[x]','aacoa[m]'},[-1 -1; 0 100]);
    model.lb(find(ismember(model.rxns,'sink_coa[c]')))=-1;
    model.ub(find(ismember(model.rxns,'sink_coa[c]')))=1;
    model.c(ismember(model.rxns,'sink_aacoa[m]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'lys-L[x] -> aacoa[m] (via Lpipecol pathway)';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 178. m8masn[r] -> nm4masn[g]
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_glc[e]'))=-1;model.ub(ismember(model.rxns,'EX_glc[e]'))=0;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=0;
    [model] = addSinkReactions(model,{'m8masn[r]','nm4masn[g]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_nm4masn[g]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'm8masn[r] -> nm4masn[g]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 179. man -> gdpmann
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'man[c]','gdpmann[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_gdpmann[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'man[c] -> gdpmann[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 180. man6p -> kdn
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'man6p[c]','kdn[c]'},[-1 -1; 0 100]);
    if ~isempty(strmatch('ACNAM9PL2',model.rxns,'exact'))
        model.c(ismember(model.rxns,'ACNAM9PL2'))=1;
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'man6p[c] -> kdn[c] - via ACNAM9PL2';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    
    %% 181. mescon[m] -> pyr[m] %changing lb has no effect
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'mescon[m]','pyr[m]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_pyr[m]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'mescon[m] -> pyr[m]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 182. met-L -> cys-L
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'met-L[c]','cys-L[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_cys-L[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'met-L[c] -> cys-L[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 183. mi145p -> inost
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'mi145p[c]','inost[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_inost[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'mi145p[c] -> inost[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 184. msa -> ala-B %changing lb has no effect
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'msa[m]','ala-B[m]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_ala-B[m]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'msa[m] -> ala-B[m]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 185. mthgxl -> 12ppd-S
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'mthgxl[c]','12ppd-S[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_12ppd-S[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'mthgxl[c] -> 12ppd-S[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 186. mthgxl -> lac-D
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'mthgxl[c]','lac-D[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_lac-D[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'mthgxl[c] -> lac-D[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 187. n2m2nmasn[l] -> man[l] + acgam[l]
    model = modelOri;
    model.c(find(model.c)) = 0;
    for i = 1 : length(RPMI_composition)
        model = changeRxnBounds(model,RPMI_composition{i},-1,'l');
    end
    [model] = addSinkReactions(model,{'n2m2nmasn[l]','man[l]','acgam[l]'},[-1 -1; 0.1 100; 0.1 100]);
    model.c(ismember(model.rxns,'sink_acgam[l]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'n2m2nmasn[l] -> man[l] + acgam[l] (with RPMI medium)';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 188. nm4masn[g] -> l2fn2m2masn[g]
    model = modelOri;
    model.c(find(model.c)) = 0;
    for i = 1 : length(RPMI_composition)
        model = changeRxnBounds(model,RPMI_composition{i},-1,'l');
    end
    [model] = addSinkReactions(model,{'nm4masn[g]','l2fn2m2masn[g]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_l2fn2m2masn[g]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'nm4masn[g] -> l2fn2m2masn[g] (with RPMI medium)';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 189. nm4masn[g] -> n2m2nmasn[g]
    model = modelOri;
    model.c(find(model.c)) = 0;
    for i = 1 : length(RPMI_composition)
        model = changeRxnBounds(model,RPMI_composition{i},-1,'l');
    end
    [model] = addSinkReactions(model,{'nm4masn[g]','n2m2nmasn[g]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_n2m2nmasn[g]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'nm4masn[g] -> n2m2nmasn[g] (with RPMI medium)';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 190. nm4masn[g] -> s2l2fn2m2masn[g]
    model = modelOri;
    model.c(find(model.c)) = 0;
    for i = 1 : length(RPMI_composition)
        model = changeRxnBounds(model,RPMI_composition{i},-1,'l');
    end
    [model] = addSinkReactions(model,{'nm4masn[g]','s2l2fn2m2masn[g]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_s2l2fn2m2masn[g]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'nm4masn[g] -> s2l2fn2m2masn[g] (with RPMI medium)';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 191. o2- -> h2o2 -> o2 + h2o
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'o2s[c]','h2o2[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_h2o2[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'o2- -> h2o2 -> o2 + h2o, 1';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    model = modelOri;
    model.c(find(model.c)) = 0;
    
    model.lb(ismember(model.rxns,'EX_h2o[e]'))=0;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=0;
    [model] = addSinkReactions(model,{'h2o2[c]','o2[c]','h2o[c]'},[-1 -1; -1 -1; 0.1 100]);
    model.c(ismember(model.rxns,'sink_h2o[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'o2- -> h2o2 -> o2 + h2o, 2';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 192. orn -> nh4 v0.05
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'orn[c]','nh4[c]'},[-1 -1; 0 100]);
    model.lb(ismember(model.rxns,'EX_nh4[e]'))=0;model.ub(ismember(model.rxns,'EX_nh4[e]'))=1000;
    model.c(ismember(model.rxns,'sink_nh4[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'orn[c] -> nh4[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 193. orn -> ptrc
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'orn[c]','ptrc[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_ptrc[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'orn[c] -> ptrc[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 194. orn -> spmd
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'orn[c]','spmd[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_spmd[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'orn[c] -> spmd[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 195. orn -> sprm
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model,rxnsInModel] = addSinkReactions(model,{'orn[c]','sprm[c]'},[-1 -1; 0 100]);
    if (rxnsInModel(2) >-1) % reaction exits already in model
        model=changeObjective(model,model.rxns(rxnsInModel(2),1));
    else
        model=changeObjective(model,'sink_sprm[c]',1);
    end
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'orn[c] -> sprm[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 196. pail_hs -> gpi_prot_hs[r]
    model = modelOri;
    model.c(find(model.c)) = 0;
    for i = 1 : length(RPMI_composition)
        model = changeRxnBounds(model,RPMI_composition{i},-1,'l');
    end
    [model] = addSinkReactions(model,{'pail_hs[c]','gpi_prot_hs[r]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_gpi_prot_hs[r]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'pail_hs[c] -> gpi_prot_hs[r] (with RMPI medium)';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 197. pail45p -> mi145p
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'pail45p_hs[c]','mi145p[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_mi145p[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'pail45p[c] -> mi145p[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 198. phe-L -> pac
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'phe-L[c]','pac[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_pac[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'phe-L[c] -> pac[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 199. phe-L -> pacald
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'phe-L[c]','pacald[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_pacald[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'phe-L[c] -> pacald[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 200. phe-L -> peamn
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'phe-L[c]','peamn[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_peamn[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'phe-L[c] -> peamn[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 201. phe-L -> phaccoa
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_glc[e]'))=-1;model.ub(ismember(model.rxns,'EX_glc[e]'))=0;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    [model] = addSinkReactions(model,{'phe-L[c]','phaccoa[c]'},[-1 -1; 0 100]);
    model.lb(find(ismember(model.rxns,'sink_phe_L[c]')))=-1;
    model.ub(find(ismember(model.rxns,'sink_phe_L[c]')))=-1;
    model.lb(find(ismember(model.rxns,'sink_coa[c]')))=-1;
    model.ub(find(ismember(model.rxns,'sink_coa[c]')))=1;
    model.c(ismember(model.rxns,'sink_phaccoa[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'phe-L[c] -> phaccoa[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 202. phe-L -> pheacgln
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'phe-L[c]','pheacgln[c]'},[-1 -1; 0 100]);
    model.lb(find(ismember(model.rxns,'sink_phe_L[c]')))=-1;
    model.ub(find(ismember(model.rxns,'sink_phe_L[c]')))=-1;
    model.c(ismember(model.rxns,'sink_pheacgln[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'phe-L[c] -> pheacgln[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 203. phe-L -> phpyr
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'phe-L[c]','phpyr[c]'},[-1 -1; 0 100]);
    model.lb(find(ismember(model.rxns,'sink_phe_L[c]')))=-1;
    model.ub(find(ismember(model.rxns,'sink_phe_L[c]')))=-1;
    model.c(ismember(model.rxns,'sink_phpyr[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'phe-L[c] -> phpyr[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 204. phe-L -> tyr-L
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'phe-L[c]','tyr-L[c]'},[-1 -1; 0 100]);
    model.lb(find(ismember(model.rxns,'sink_phe_L[c]')))=-1;
    model.ub(find(ismember(model.rxns,'sink_phe_L[c]')))=-1;
    model.c(ismember(model.rxns,'sink_tyr-L[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'phe-L[c] -> tyr-L[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 205. pheme -> bilirub %changing lb has no effect
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_glc[e]'))=-1;model.ub(ismember(model.rxns,'EX_glc[e]'))=0;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    [model] = addSinkReactions(model,{'pheme[c]','bilirub[c]'},[-1 -1; 0 100]);
    model.lb(find(ismember(model.rxns,'sink_pheme[c]')))=-1;
    model.ub(find(ismember(model.rxns,'sink_pheme[c]')))=-1;
    model.c(ismember(model.rxns,'sink_bilirub[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'pheme[c] -> bilirub[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 205. phytcoa[x] -> dmnoncoa[m]
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'phytcoa[x]','dmnoncoa[m]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_dmnoncoa[m]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'phytcoa[x] -> dmnoncoa[m]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 206. pmtcoa[c] -> crmp_hs[c]
    model = modelOri;
    for i = 1 : length(RPMI_composition)
        model = changeRxnBounds(model,RPMI_composition{i},-1,'l');
    end
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'pmtcoa[c]','crmp_hs[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_crmp_hs[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'pmtcoa[c] -> crmp_hs[c] (with RPMI medium)';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 207. pmtcoa[c] -> sphmyln_hs[c]
    model = modelOri;
    for i = 1 : length(RPMI_composition)
        model = changeRxnBounds(model,RPMI_composition{i},-1,'l');
    end
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'pmtcoa[c]','sphmyln_hs[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_sphmyln_hs[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'pmtcoa[c] -> sphmyln_hs[c] (with RPMI medium)';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 208. ppcoa[m] -> succoa[m]
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'ppcoa[m]','succoa[m]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_succoa[m]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'ppcoa[m] -> succoa[m]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 209. pro-L -> glu-L
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'pro-L[c]','glu-L[c]'},[-1 -1; 0 100]);
    model.lb(find(ismember(model.rxns,'sink_pro_L[c]')))=-1;
    model.ub(find(ismember(model.rxns,'sink_pro_L[c]')))=-1;
    model.c(ismember(model.rxns,'sink_glu-L[c]'))=1;
    model.c(ismember(model.rxns,'sink_glu_L[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        if isempty(eFBA.f)
            eFBA.f=0;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'pro-L[c] -> glu-L[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 210. ptrc -> ala-B
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'ptrc[c]','ala-B[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_ala-B[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'ptrc[c] -> ala-B[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 211. ptrc -> spmd
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'ptrc[c]','spmd[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_spmd[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'ptrc[c] -> spmd[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 212. pyr -> fad[m] + h[m] %changing lb has no effect
    model = modelOri;
    model.c(find(model.c)) = 0;
    for i = 1 : length(RPMI_composition)
        model = changeRxnBounds(model,RPMI_composition{i},-1,'l');
    end
    [model] = addSinkReactions(model,{'pyr[c]','fadh2[m]','fad[m]','h[m]'},[-1 -1;-1 0;  0.1 100; 0.1 100]);
    model.c(ismember(model.rxns,'sink_fad[m]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'pyr -> fad[m] + h[m] (with RPMI medium)';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 213. pyr -> lac-D
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'pyr[c]','lac-D[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_lac-D[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'pyr[c] -> lac-D[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 214. pyr -> nad[m] + h[m] %changing lb has no effect
    model = modelOri;
    model.c(find(model.c)) = 0;
    for i = 1 : length(RPMI_composition)
        model = changeRxnBounds(model,RPMI_composition{i},-1,'l');
    end
    [model] = addSinkReactions(model,{'pyr[c]','nad[m]','h[m]'},[-1 -1; 0.1 100; 0.1 100]);
    model.c(ismember(model.rxns,'sink_nad[m]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'pyr -> nad[m] + h[m] (with RPMI medium)';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 215. pyr[c] -> accoa[m] + co2[c] + nadh[m] %changing lb has no effect
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'pyr[c]','accoa[m]','nadh[m]','co2[c]'},[-1 -1; 0.1 100; 0.1 100; 0.1 100]);
    model.lb(find(ismember(model.rxns,'sink_coa[c]')))=-1;
    model.ub(find(ismember(model.rxns,'sink_coa[c]')))=1;
    model.lb(find(ismember(model.rxns,'sink_nad[c]')))=-1;
    model.ub(find(ismember(model.rxns,'sink_nad[c]')))=1;
    model.c(ismember(model.rxns,'sink_accoa[m]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'pyr[c] -> accoa[m] + co2[c] + nadh[m]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 216. pyr<>ala-L
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'pyr[c]','ala-L[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_ala-L[c]'))=1;
    model.c(ismember(model.rxns,'sink_ala_L[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'pyr[c] -> ala-L[c], 1';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'ala-L[c]','pyr[c]'},[-1 -1; 0 100]);
    model.lb(find(ismember(model.rxns,'sink_ala_L[c]')))=-1;
    model.ub(find(ismember(model.rxns,'sink_ala_L[c]')))=-1;
    model.c(ismember(model.rxns,'sink_pyr[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'pyr[c] -> ala-L[c], 2';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% R_group
    %% 217. s2l2fn2m2masn[l] -> man[l] + acgam[l]
    model = modelOri;
    model.c(find(model.c)) = 0;
    for i = 1 : length(RPMI_composition)
        model = changeRxnBounds(model,RPMI_composition{i},-1,'l');
    end
    [model] = addSinkReactions(model,{'s2l2fn2m2masn[l]','man[l]','acgam[l]'},[-1 -1; 0.1 100; 0.1 100]);
    model.c(ismember(model.rxns,'sink_man[l]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 's2l2fn2m2masn[l] -> man[l] + acgam[l] (with RPMI medium)';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 218. Ser/Thr[g] + udpacgal[g] -> core2[g] %changing lb has no effect
    model = modelOri;
    model.c(find(model.c)) = 0;
    for i = 1 : length(RPMI_composition)
        model = changeRxnBounds(model,RPMI_composition{i},-1,'l');
    end
    model = changeRxnBounds(model,'GALNTg',0.1,'l');
    [model] = addSinkReactions(model,{'Ser/Thr[g]';'udpacgal[g]';'core2[g]'},[-1 -1; -1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_core2[g]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Ser/Thr[g] + udpacgal[g] -> core2[g] - via GALNTg and DM_core4[g] (with RPMI medium)';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 219. Ser/Thr[g] + udpacgal[g] -> core4[g] %changing lb has no effect
    model = modelOri;
    model.c(find(model.c)) = 0;
    for i = 1 : length(RPMI_composition)
        model = changeRxnBounds(model,RPMI_composition{i},-1,'l');
    end
    model = changeRxnBounds(model,'GALNTg',0.1,'l');
    [model] = addSinkReactions(model,{'Ser/Thr[g]';'udpacgal[g]';'core4[g]'},[-1 -1; -1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_core4[g]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Ser/Thr[g] + udpacgal[g] -> core4[g] - via GALNTg and DM_core4[g] (with RPMI medium)';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 220. Ser/Thr[g] + udpacgal[g] -> Tn_antigen[g] % dsTn_antigen does not exists
    % - I used Tn_antigen instead %changing lb has no effect
    model = modelOri;
    model.c(find(model.c)) = 0;
    for i = 1 : length(RPMI_composition)
        model = changeRxnBounds(model,RPMI_composition{i},-1,'l');
    end
    [model] = addSinkReactions(model,{'Ser/Thr[g]','udpacgal[g]','Tn_antigen[g]'},[-1 -1;-1 -1; 0 100]);
    % model.c(ismember(model.rxns,'sink_Tn_antigen[g]'))=1;
    if ~isempty(strmatch('GALNTg',model.rxns,'exact'))
        model = changeObjective(model, 'GALNTg',1);
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Ser/Thr[g] + udpacgal[g] -> Tn_antigen[g] - via GALNTg (with RPMI medium)';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    
    %% 221. Ser/Thr[g] + udpacgal[g] -> sTn_antigen[g] %changing lb has no effect
    model = modelOri;
    model.c(find(model.c)) = 0;
    for i = 1 : length(RPMI_composition)
        model = changeRxnBounds(model,RPMI_composition{i},-1,'l');
    end
    if ~isempty(strmatch('GALNTg',model.rxns,'exact'))
        model = changeRxnBounds(model,'GALNTg',0.1,'l');
        [model] = addSinkReactions(model,{'Ser/Thr[g]','udpacgal[g]','sTn_antigen[g]'},[-1 -1;-1 -1; 0 100]);
        if (rxnsInModel(1) >-1) % reaction exits already in model
            model=changeObjective(model,model.rxns(rxnsInModel(1),1));
        else
            model=changeObjective(model,'sink_sTn_antigen[g]',1);
        end
        if find(model.c)>0
            try
                eFBA=entropicFluxBalanceAnalysis(model,param);
                eFBA.f=sum(eFBA.v(find(model.c>0)));
            catch
                eFBA.f=NaN;
            end
            TestSolution(k,1) = eFBA.f;
        else
            TestSolution(k,1) = NaN;
        end
        TestSolutionName{k,1} = 'Ser/Thr[g] + udpacgal[g] -> sTn_antigen[g] - via GALNTg and DM_sTn_antigen[g] (with RPMI medium)';
        if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    else
        TestSolution(k,1) = NaN;
        TestSolutionName{k,1} = 'Ser/Thr[g] + udpacgal[g] -> sTn_antigen[g] - via GALNTg and DM_sTn_antigen[g]';
        if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    end
    %% 222. Ser-Gly/Ala-X-Gly[er] -> cs_pre[g]
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_glc[e]'))=-1;model.ub(ismember(model.rxns,'EX_glc[e]'))=0;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    [model] = addSinkReactions(model,{'Ser-Gly/Ala-X-Gly[r]','cs_pre[g]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_cs_pre[g]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Ser-Gly/Ala-X-Gly[r] -> cs_pre[g]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 223. Ser-Gly/Ala-X-Gly[er] -> cspg_a[g]
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_glc[e]'))=-1;model.ub(ismember(model.rxns,'EX_glc[e]'))=0;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    [model] = addSinkReactions(model,{'Ser-Gly/Ala-X-Gly[r]','cspg_a[g]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_cspg_a[g]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Ser-Gly/Ala-X-Gly[r] -> cspg_a[g]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 224. Ser-Gly/Ala-X-Gly[er] -> cspg_c[g]
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_glc[e]'))=-1;model.ub(ismember(model.rxns,'EX_glc[e]'))=0;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    [model] = addSinkReactions(model,{'Ser-Gly/Ala-X-Gly[r]','cspg_c[g]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_cspg_c[g]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Ser-Gly/Ala-X-Gly[r] -> cspg_c[g]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 225. Ser-Gly/Ala-X-Gly[er] -> cspg_d[g]
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_glc[e]'))=-1;model.ub(ismember(model.rxns,'EX_glc[e]'))=0;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    [model] = addSinkReactions(model,{'Ser-Gly/Ala-X-Gly[r]','cspg_d[g]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_cspg_d[g]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Ser-Gly/Ala-X-Gly[r] -> cspg_d[g]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 226. Ser-Gly/Ala-X-Gly[er] -> cspg_e[g]
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_glc[e]'))=-1;model.ub(ismember(model.rxns,'EX_glc[e]'))=0;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    [model] = addSinkReactions(model,{'Ser-Gly/Ala-X-Gly[r]','cspg_e[g]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_cspg_e[g]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Ser-Gly/Ala-X-Gly[r] -> cspg_e[g]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 227. Ser-Gly/Ala-X-Gly[er] -> hspg[g]
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_glc[e]'))=-1;model.ub(ismember(model.rxns,'EX_glc[e]'))=0;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    [model] = addSinkReactions(model,{'Ser-Gly/Ala-X-Gly[r]','hspg[g]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_hspg[g]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Ser-Gly/Ala-X-Gly[r] -> hspg[g]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 228. Ser-Gly/Ala-X-Ser[er] -> cspg_b[g]
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_glc[e]'))=-1;model.ub(ismember(model.rxns,'EX_glc[e]'))=0;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    [model] = addSinkReactions(model,{'Ser-Gly/Ala-X-Gly[r]','cspg_b[g]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_cspg_b[g]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Ser-Gly/Ala-X-Gly[r] -> cspg_b[g]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 229. ser-L -> cys-L
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_glc[e]'))=-1;model.ub(ismember(model.rxns,'EX_glc[e]'))=0;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    [model] = addSinkReactions(model,{'ser-L[c]','cys-L[c]'},[-1 -1; 0 100]);
    model.lb(find(ismember(model.rxns,'sink_ser_L[c]')))=-1;
    model.ub(find(ismember(model.rxns,'sink_ser_L[c]')))=-1;
    model.c(ismember(model.rxns,'sink_cys-L[c]'))=1;
    model.c(ismember(model.rxns,'sink_cys_L[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'ser-L[c] -> cys-L[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 230. so4 -> PAPS
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_glc[e]'))=-1;model.ub(ismember(model.rxns,'EX_glc[e]'))=0;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    [model] = addSinkReactions(model,{'so4[c]','paps[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_paps[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'so4[c] -> paps[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 231. spmd -> sprm
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model,rxnsInModel] = addSinkReactions(model,{'spmd[c]','sprm[c]'},[-1 -1; 0 100]);
    %  model.c(ismember(model.rxns,'sink_sprm[c]'))=1;
    if (rxnsInModel(2) >-1) % reaction exits already in model
        model=changeObjective(model,model.rxns(rxnsInModel(2),1));
    else
        model=changeObjective(model,'sink_sprm[c]',1);
    end
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'spmd[c] -> sprm[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 232. srtn -> f5hoxkyn
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_glc[e]'))=-1;model.ub(ismember(model.rxns,'EX_glc[e]'))=0;
    [model] = addSinkReactions(model,{'srtn[c]','f5hoxkyn[c]'},[-1 -1; 0 100]);
    model.lb(find(ismember(model.rxns,'DM_srtn[c]')))=-1;
    model.ub(find(ismember(model.rxns,'DM_srtn[c]')))=-1;
    model.c(ismember(model.rxns,'sink_f5hoxkyn[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'srtn[c] -> f5hoxkyn[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 233. srtn -> fna5moxam
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_glc[e]'))=-1;model.ub(ismember(model.rxns,'EX_glc[e]'))=0;
    [model] = addSinkReactions(model,{'srtn[c]','fna5moxam[c]'},[-1 -1; 0 100]);
    model.lb(find(ismember(model.rxns,'DM_srtn[c]')))=-1;
    model.ub(find(ismember(model.rxns,'DM_srtn[c]')))=-1;
    model.c(ismember(model.rxns,'sink_fna5moxam[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'srtn[c] -> fna5moxam[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 234. srtn -> nmthsrtn
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'srtn[c]','nmthsrtn[c]'},[-1 -1; 0 100]);
    model.lb(find(ismember(model.rxns,'DM_srtn[c]')))=-1;
    model.ub(find(ismember(model.rxns,'DM_srtn[c]')))=-1;
    model.c(find(ismember(model.rxns,'sink_nmthsrtn[c]')))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'srtn[c] -> nmthsrtn[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 235. strch1[e] -> glc[e]
    model = modelOri;
    model.c(find(model.c)) = 0;
    model = changeRxnBounds(model,'EX_strch1[e]',-1,'l');
    model = changeRxnBounds(model,'EX_strch1[e]',-1,'u');
    model = changeRxnBounds(model,'EX_glc[e]',0,'l');
    model = changeRxnBounds(model,'EX_glc[e]',1000,'u');
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=0;
    if ~isempty(strmatch('AMY1e',model.rxns,'exact'))
        model.c(ismember(model.rxns,'AMY1e'))=1;
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'strch1[e] -> glc-D[e] via AMY1e';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    
    %% 236. succoa[m] -> oaa[m]
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'succoa[m]','oaa[m]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_oaa[m]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'succoa[m] -> oaa[m]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 237. taur[x] -> tchola[x]
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_glc[e]'))=-1;model.ub(ismember(model.rxns,'EX_glc[e]'))=0;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    [model] = addSinkReactions(model,{'taur[x]','tchola[x]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_tchola[x]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'taur[x] -> tchola[x]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 238. thcholstoic[x] -> gchola[x]
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_glc[e]'))=-1;model.ub(ismember(model.rxns,'EX_glc[e]'))=0;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    [model] = addSinkReactions(model,{'thcholstoic[x]','gchola[x]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_gchola[x]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'thcholstoic[x] -> gchola[x]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 239. thcholstoic[x] -> tchola[x]
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_glc[e]'))=-1;model.ub(ismember(model.rxns,'EX_glc[e]'))=0;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    [model] = addSinkReactions(model,{'thcholstoic[x]','tchola[x]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_tchola[x]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'thcholstoic[x] -> tchola[x]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 240. thr-L -> ppcoa
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    [model] = addSinkReactions(model,{'trp-L[c]','ppcoa[c]'},[-1 -1; 0 100]);
    model.lb(find(ismember(model.rxns,'sink_trp_L[c]')))=-1;
    model.ub(find(ismember(model.rxns,'sink_trp_L[c]')))=-1;
    model.lb(find(ismember(model.rxns,'sink_coa[c]')))=-1;
    model.ub(find(ismember(model.rxns,'sink_coa[c]')))=1;
    model.c(ismember(model.rxns,'sink_ppcoa[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'trp-L[c] -> ppcoa[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 241. trp-L -> accoa
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    [model] = addSinkReactions(model,{'trp-L[c]','accoa[c]'},[-1 -1; 0 100]);
    model.lb(find(ismember(model.rxns,'sink_trp_L[c]')))=-1;
    model.ub(find(ismember(model.rxns,'sink_trp_L[c]')))=-1;
    model.lb(find(ismember(model.rxns,'sink_coa[c]')))=-1;
    model.ub(find(ismember(model.rxns,'sink_coa[c]')))=1;
    model.c(ismember(model.rxns,'sink_accoa[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'trp-L[c] -> accoa[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 242. trp-L -> anth
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model,rxnsInModel] = addSinkReactions(model,{'trp-L[c]','anth[c]'},[-1 -1; 0 100]);
    model.lb(find(ismember(model.rxns,'sink_trp_L[c]')))=-1;
    model.ub(find(ismember(model.rxns,'sink_trp_L[c]')))=-1;
    if (rxnsInModel(2) >-1) % reaction exits already in model
        model=changeObjective(model,model.rxns(rxnsInModel(2),1));
    else
        model=changeObjective(model,'sink_anth[c]',1);
    end
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'trp-L[c] -> anth[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 243. trp-L -> id3acald
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'trp-L[c]','id3acald[c]'},[-1 -1; 0 100]);
    model.lb(find(ismember(model.rxns,'sink_trp_L[c]')))=-1;
    model.ub(find(ismember(model.rxns,'sink_trp_L[c]')))=-1;
    model.c(ismember(model.rxns,'sink_id3acald[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'trp-L[c] -> id3acald[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 244. trp-L -> kynate
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'trp-L[c]','kynate[c]'},[-1 -1; 0 100]);
    model.lb(find(ismember(model.rxns,'sink_trp_L[c]')))=-1;
    model.ub(find(ismember(model.rxns,'sink_trp_L[c]')))=-1;
    model.c(ismember(model.rxns,'sink_kynate[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'trp-L[c] -> kynate[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 245. trp-L -> melatn
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'trp-L[c]','melatn[c]'},[-1 -1; 0 100]);
    model.lb(find(ismember(model.rxns,'sink_trp_L[c]')))=-1;
    model.ub(find(ismember(model.rxns,'sink_trp_L[c]')))=-1;
    model.c(ismember(model.rxns,'sink_melatn[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'trp-L[c] -> melatn[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 246. trp-L -> melatn
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'trp-L[c]','Lfmkynr[c]'},[-1 -1; 0 100]);
    model.lb(find(ismember(model.rxns,'sink_trp_L[c]')))=-1;
    model.ub(find(ismember(model.rxns,'sink_trp_L[c]')))=-1;
    model.c(ismember(model.rxns,'sink_Lfmkynr[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'trp-L[c] -> Lfmkynr[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 247. trp-L -> melatn
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'trp-L[c]','Lkynr[c]'},[-1 -1; 0 100]);
    model.lb(find(ismember(model.rxns,'sink_trp_L[c]')))=-1;
    model.ub(find(ismember(model.rxns,'sink_trp_L[c]')))=-1;
    model.c(ismember(model.rxns,'sink_Lkynr[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'trp-L[c] -> Lkynr[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 248. trp-L -> melatn
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'trp-L[c]','nformanth[c]'},[-1 -1; 0 100]);
    model.lb(find(ismember(model.rxns,'sink_trp_L[c]')))=-1;
    model.ub(find(ismember(model.rxns,'sink_trp_L[c]')))=-1;
    model.c(ismember(model.rxns,'sink_nformanth[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'trp-L[c] -> nformanth[c]';
    model.lb(find(ismember(model.rxns,'sink_trp_L[c]')))=-1;
    model.ub(find(ismember(model.rxns,'sink_trp_L[c]')))=-1;
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 249. srtn[c] -> 5moxact[c]
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_glc[e]'))=-1;model.ub(ismember(model.rxns,'EX_glc[e]'))=0;
    [model] = addSinkReactions(model,{'srtn[c]','5moxact[c]'},[-1 -1; 0 100]);
    model.lb(find(ismember(model.rxns,'DM_srtn[c]')))=-1;
    model.ub(find(ismember(model.rxns,'DM_srtn[c]')))=-1;
    model.c(ismember(model.rxns,'sink_5moxact[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'srtn[c] -> 5moxact[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 250. srtn[c] -> 6hoxmelatn[c]
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_glc[e]'))=-1;model.ub(ismember(model.rxns,'EX_glc[e]'))=0;
    [model] = addSinkReactions(model,{'srtn[c]','6hoxmelatn[c]'},[-1 -1; 0 100]);
    model.lb(find(ismember(model.rxns,'DM_srtn[c]')))=-1;
    model.ub(find(ismember(model.rxns,'DM_srtn[c]')))=-1;
    model.c(ismember(model.rxns,'sink_6hoxmelatn[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'srtn[c] -> 6hoxmelatn[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 251. trp-L -> quln
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'trp-L[c]','quln[c]'},[-1 -1; 0 100]);
    model.lb(find(ismember(model.rxns,'sink_trp_L[c]')))=-1;
    model.ub(find(ismember(model.rxns,'sink_trp_L[c]')))=-1;
    model.c(ismember(model.rxns,'sink_quln[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'trp-L[c] -> quln[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 252. trp-L -> srtn
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'trp-L[c]','srtn[c]'},[-1 -1; 0 100]);
    model.lb(find(ismember(model.rxns,'sink_trp_L[c]')))=-1;
    model.ub(find(ismember(model.rxns,'sink_trp_L[c]')))=-1;
    model.c(ismember(model.rxns,'DM_srtn[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'trp-L[c] -> srtn[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 253. Tyr-ggn -> glygn2
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    [model] = addSinkReactions(model,{'Tyr-ggn[c]','glygn2[c]'},[-1 -1; 0 100]);
    model.lb(find(ismember(model.rxns,'sink_Tyr-ggn[c]')))=-1;
    model.ub(find(ismember(model.rxns,'sink_Tyr-ggn[c]')))=-1;
    model.c(ismember(model.rxns,'sink_glygn2[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Tyr-ggn[c] -> glygn2[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 254. tyr-L -> 34hpp
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'tyr-L[c]','34hpp[c]'},[-1 -1; 0 100]);
    model.lb(find(ismember(model.rxns,'sink_tyr_L[c]')))=-1;
    model.ub(find(ismember(model.rxns,'sink_tyr_L[c]')))=-1;
    model.c(ismember(model.rxns,'sink_34hpp[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'tyr-L[c] -> 34hpp[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 255. tyr-L -> 4hphac
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'tyr-L[c]','4hphac[c]'},[-1 -1; 0 100]);
    model.lb(find(ismember(model.rxns,'sink_tyr_L[c]')))=-1;
    model.ub(find(ismember(model.rxns,'sink_tyr_L[c]')))=-1;
    model.c(ismember(model.rxns,'sink_4hphac[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'tyr-L[c] -> 4hphac[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 256. tyr-L -> adrnl
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'tyr-L[c]','adrnl[c]'},[-1 -1; 0 100]);
    model.lb(find(ismember(model.rxns,'sink_tyr_L[c]')))=-1;
    model.ub(find(ismember(model.rxns,'sink_tyr_L[c]')))=-1;
    model.c(ismember(model.rxns,'sink_adrnl[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'tyr-L[c] -> adrnl[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 257. tyr-L -> dopa
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'tyr-L[c]','dopa[c]'},[-1 -1; 0 100]);
    model.lb(find(ismember(model.rxns,'sink_tyr_L[c]')))=-1;
    model.ub(find(ismember(model.rxns,'sink_tyr_L[c]')))=-1;
    model.c(ismember(model.rxns,'sink_dopa[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'tyr-L[c] -> dopa[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 258. tyr-L -> fum + acac
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'tyr-L[c]','fum[c]','acac[c]'},[-1 -1; 0.1 100; 0.1 100]);
    model.lb(find(ismember(model.rxns,'sink_tyr_L[c]')))=-1;
    model.ub(find(ismember(model.rxns,'sink_tyr_L[c]')))=-1;
    model.c(ismember(model.rxns,'sink_fum[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'tyr-L[c] -> fum[c] + acac[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 259. tyr-L -> melanin
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model,rxnsInModel] = addSinkReactions(model,{'tyr-L[c]','melanin[c]'},[-1 -1; 0 100]);
    model.lb(find(ismember(model.rxns,'sink_tyr_L[c]')))=-1;
    model.ub(find(ismember(model.rxns,'sink_tyr_L[c]')))=-1;
    if (rxnsInModel(2) >-1) % reaction exits already in model
        model=changeObjective(model,model.rxns(rxnsInModel(2),1));
    else
        model=changeObjective(model,'sink_melanin[c]',1);
    end
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'tyr-L[c] -> melanin[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 260. tyr-L -> nrpphr
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'tyr-L[c]','nrpphr[c]'},[-1 -1; 0 100]);
    model.lb(find(ismember(model.rxns,'sink_tyr_L[c]')))=-1;
    model.ub(find(ismember(model.rxns,'sink_tyr_L[c]')))=-1;
    model.c(ismember(model.rxns,'sink_nrpphr[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'tyr-L[c] -> nrpphr[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 261. uacgam + udpglcur -> ha[e] %changing lb has no effect
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_glc[e]'))=-1;model.ub(ismember(model.rxns,'EX_glc[e]'))=0;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    [model] = addSinkReactions(model,{'uacgam[c]','udpglcur[c]','ha[e]'},[-1 -1; -1 -1;0 100]);
    %model.c(ismember(model.rxns,'sink_ha[c]'))=1;
    if ~isempty(strmatch('HAS2',model.rxns,'exact'))
        model=changeObjective(model,'HAS2',1);
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} =  'uacgamv[c] + udpglcur[c] -> ha[e] - via HAS2';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA   
    %% 262. uacgam -> m8masn[r]
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    [model] = addSinkReactions(model,{'uacgam[c]','m8masn[r]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_m8masn[r]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'uacgam[c] -> m8masn[r]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 263. udpglcur -> xu5p-D
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'udpglcur[c]','xu5p-D[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_xu5p-D[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'udpglcur[c] -> xu5p-D[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 264. ura -> ala-B
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'ura[c]','ala-B[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_ala-B[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'ura[c] -> ala-B[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 265. val-L -> 3aib
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'val-L[c]','3aib[c]'},[-1 -1; 0 100]);
    model.lb(find(ismember(model.rxns,'sink_val_L[c]')))=-1;
    model.ub(find(ismember(model.rxns,'sink_val_L[c]')))=-1;
    model.c(ismember(model.rxns,'sink_3aib[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'val-L[c] -> 3aib[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 266. val-L -> succoa
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    [model] = addSinkReactions(model,{'val-L[c]','succoa[m]'},[-1 -1; 0 100]);
    model.lb(find(ismember(model.rxns,'sink_val_L[c]')))=-1;
    model.ub(find(ismember(model.rxns,'sink_val_L[c]')))=-1;
    model.lb(find(ismember(model.rxns,'sink_coa[c]')))=-1;
    model.ub(find(ismember(model.rxns,'sink_coa[c]')))=1;
    model.c(ismember(model.rxns,'sink_succoa[m]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'val-L[c] -> succoa[m]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 267. xoltriol[m] -> thcholstoic[m]
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'xoltriol[m]','thcholstoic[m]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_thcholstoic[m]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'xoltriol[m] -> thcholstoic[m]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 268. xylu-D -> glyclt
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'xylu-D[c]','glyclt[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_glyclt[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'xylu-D[c] -> glyclt[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
end

%% metabolic tasks based on Enterocyte model - without original ('u')
% compartment I deleted the last argument in changeObjective from here
% onwards (SS)
if strcmp(test,'IEC') || strcmp(test,'all')|| strcmp(test,'Harvey')
    %% 269. glucose to lactate conversion
    model=modelOri;
    model=changeRxnBounds(model,'EX_glc[e]',-1,'b');
    if ~isempty(strmatch('EX_lac-L[e]',model.rxns,'exact'))
        model=changeObjective(model,'EX_lac-L[e]',1);
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'glucose to lactate conversion';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    
    %%  270. glutamine to glucose conversion
    model=modelOri;
    model=changeRxnBounds(model,'EX_gln-L[e]',-1,'b');
    model=changeRxnBounds(model,'EX_glc[e]',0,'b');
    model=changeRxnBounds(model,'EX_malt[e]',0,'b');
    model=changeRxnBounds(model,'EX_strch1[e]',0,'b');
    model=changeRxnBounds(model,'EX_strch2[e]',0,'b');
    model=changeRxnBounds(model,'EX_sucr[e]',0,'b');
    
    if ~isempty(strmatch('GLUNm',model.rxns,'exact'))
        model=changeObjective(model,'GLUNm',1);
        %eFBA=optimizeCbModel(model,'min');
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'glutamine to glucose conversion - GLUNm';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    
    %% 271. glutamine to glucose conversion - ASPTAm
    if ~isempty(strmatch('ASPTAm',model.rxns,'exact'))
        model=changeObjective(model,'ASPTAm',1);
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'glutamine to glucose conversion - ASPTAm';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    
    %% 272. 'glutamine to glucose conversion - FUM'
    if ~isempty(strmatch('FUM',model.rxns,'exact'))
        model=changeObjective(model,'FUM',1);
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'glutamine to glucose conversion - FUM';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 273. glutamine to glucose conversion - MDH
    if ~isempty(strmatch('MDH',model.rxns,'exact'))
        model=changeObjective(model,'MDH',1);
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'glutamine to glucose conversion - MDH';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 274. glutamine to glucose conversion - G6PPer
    if ~isempty(strmatch('G6PPer',model.rxns,'exact'))
        model=changeObjective(model,'G6PPer',1);
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'glutamine to glucose conversion - G6PPer';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 274. glutamine to proline conversion
    model=modelOri;
    model=changeRxnBounds(model,'EX_gln-L[e]',-1,'b');
    model=changeRxnBounds(model,'EX_glc[e]',0,'b');
    model=changeRxnBounds(model,'EX_pro-L[e]',0,'b');
    
    if ~isempty(strmatch('P5CRm',model.rxns,'exact'))
        model=changeObjective(model,'P5CRm',1);
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'glutamine to proline conversion - P5CRm';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 275. glutamine to proline conversion - P5CRxm
    model=modelOri;
    model=changeRxnBounds(model,'EX_gln-L[e]',-1,'b');
    model=changeRxnBounds(model,'EX_glc[e]',0,'b');
    model=changeRxnBounds(model,'EX_pro-L[e]',0,'b');
    
    if ~isempty(strmatch('P5CRxm',model.rxns,'exact'))
        model=changeObjective(model,'P5CRxm',1);
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'glutamine to proline conversion - P5CRxm';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA 
    %% 276. glutamine to ornithine conversion
    model=modelOri;
    
    model=changeRxnBounds(model,'EX_gln-L[e]',-1,'l');
    if ~isempty(strmatch('ORNTArm',model.rxns,'exact'))
        model=changeObjective(model,'ORNTArm',1);
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'glutamine to ornithine conversion - ORNTArm';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    
    %% 277. glutamine to citrulline converion
    model=modelOri;
    model=changeRxnBounds(model,'EX_gln-L[e]',-1,'b');
    if ~isempty(strmatch('OCBTm',model.rxns,'exact'))
        model=changeObjective(model,'OCBTm',1);
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'glutamine to citrulline converion - OCBTm';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    
    %% 278. glutamine to lactate
    model=modelOri;
    model=changeRxnBounds(model,'EX_gln-L[e]',-1,'b');
    if ~isempty(strmatch('LDH_L',model.rxns,'exact'))
        model=changeObjective(model,'LDH_L',1);
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'glutamine to lactate - LDH_L';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    
    %% 279. glutamine to aspartate
    model=modelOri;
    model=changeRxnBounds(model,'EX_gln-L[e]',-1,'b');
    if ~isempty(strmatch('ASPTA',model.rxns,'exact'))
        model=changeObjective(model,'ASPTA',1);
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'glutamine to aspartate - ASPTA';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    
    %% 280. glutamine to co2
    model=modelOri;
    model=changeRxnBounds(model,'EX_gln-L[e]',-1,'b');
    if ~isempty(strmatch('AKGDm',model.rxns,'exact'))
        model=changeObjective(model,'AKGDm',1);
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'glutamine to co2 - AKGDm';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 281. glutamine to ammonia
    model=modelOri;
    model=changeRxnBounds(model,'EX_gln-L[e]',-1,'b');
    if ~isempty(strmatch('GLUNm',model.rxns,'exact'))
        model=changeObjective(model,'GLUNm',1);
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'glutamine to ammonia - GLUNm';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    
    %% 282. putriscine to methionine (depends on oxygen uptake);
    model=modelOri;
    model=changeRxnBounds(model,'EX_ptrc[e]',-1,'b');
    model=changeRxnBounds(model,'EX_o2[e]',-1,'l');
    if ~isempty(strmatch('UNK2',model.rxns,'exact'))
        model=changeObjective(model,'UNK2',1);
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'putriscine to methionine (depends on oxygen uptake) - UNK2';
    if ~isnan(TestSolution(k,1)); if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ; end ;k = k +1;clear eFBA
    
    %%  283. secretion of alanine
    model=modelOri;
    
    model.lb(ismember(model.rxns,'EX_glc[e]'))=-1;model.ub(ismember(model.rxns,'EX_glc[e]'))=0;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    if ~isempty(strmatch('EX_ala_L[e]',model.rxns,'exact'))
        model=changeObjective(model,'EX_ala_L[e]',1);
        %eFBA=optimizeCbModel(model,'min');
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'secretion of alanine';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    
    %%  284. secretion of lactate
    model=modelOri;
    model.lb(ismember(model.rxns,'EX_glc[e]'))=-1;model.ub(ismember(model.rxns,'EX_glc[e]'))=0;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    if ~isempty(strmatch('EX_lac-L[e]',model.rxns,'exact'))
        model=changeObjective(model,'EX_lac-L[e]');
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'secretion of lactate';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    
    %% 285. synthesis of arginine from glutamine
    model=modelOri;
    model=changeRxnBounds(model,'EX_gln-L[e]',-1,'b');
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    if ~isempty(strmatch('ARGSL',model.rxns,'exact'))
        model=changeObjective(model,'ARGSL',1);
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'synthesis of arginine from glutamine - ARGSL';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA 
    %% 286. synthesis of proline from glutamine
    model=modelOri;
    model=changeRxnBounds(model,'EX_gln-L[e]',-1,'b');
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    if ~isempty(strmatch('P5CR',model.rxns,'exact'))
        model=changeObjective(model,'P5CR',1);
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'synthesis of proline from glutamine - P5CR';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 287. synthesis of proline from glutamine
    model=modelOri;
    model=changeRxnBounds(model,'EX_gln-L[e]',-1,'b');
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    if ~isempty(strmatch('P5CRm',model.rxns,'exact'))
        model=changeObjective(model,'P5CRm',1);
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'synthesis of proline from glutamine - P5CRm';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    
    %% 288. synthesis of proline from glutamine
    model=modelOri;
    model=changeRxnBounds(model,'EX_gln-L[e]',-1,'b');
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    if ~isempty(strmatch('P5CRxm',model.rxns,'exact'))
        model=changeObjective(model,'P5CRxm',1);
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'synthesis of proline from glutamine - P5CRxm';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    
    %% 289. synthesis of alanine from glutamine
    model=modelOri;
    model=changeRxnBounds(model,'EX_o2[e]',-1,'l');
    model=changeRxnBounds(model,'EX_gln-L[e]',-1,'b');
    if ~isempty(strmatch('ALATA_L',model.rxns,'exact'))
        model=changeObjective(model,'ALATA_L',1);
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'synthesis of alanine from glutamine - ALATA_L';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    
    %% 290. basolateral secretion of proline
    model=modelOri;
    model.lb(ismember(model.rxns,'EX_glc[e]'))=-1;model.ub(ismember(model.rxns,'EX_glc[e]'))=0;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    if ~isempty(strmatch('EX_pro-L[e]',model.rxns,'exact'))
        model=changeObjective(model,'EX_pro-L[e]',1);
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'secretion of proline';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    
    %% 291. basolateral secretion of arginine
    model=modelOri;
    model.lb(ismember(model.rxns,'EX_glc[e]'))=-1;model.ub(ismember(model.rxns,'EX_glc[e]'))=0;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    if ~isempty(strmatch('EX_arg-L[e]',model.rxns,'exact'))
        model=changeObjective(model,'EX_arg-L[e]',1);
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'secretion of arginine';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    
    %% 292. basolateral secretion of ornithine
    model=modelOri;
    model.lb(ismember(model.rxns,'EX_glc[e]'))=-1;model.ub(ismember(model.rxns,'EX_glc[e]'))=0;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    if ~isempty(strmatch('EX_orn[e]',model.rxns,'exact'))
        model=changeObjective(model,'EX_orn[e]',1);
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'secretion of ornithine';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 293. synthesis of spermine from ornithine
    model=modelOri;
    model.lb(ismember(model.rxns,'EX_orn[e]'))=-1;model.ub(ismember(model.rxns,'EX_orn[e]'))=-1;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    if ~isempty(strmatch('SPRMS',model.rxns,'exact'))
        model=changeObjective(model,'SPRMS',1);
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'synthesis of spermine from ornithine - SPRMS';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 294. synthesis of spermidine from ornithine
    model=modelOri;
    model.lb(ismember(model.rxns,'EX_orn[e]'))=-1;model.ub(ismember(model.rxns,'EX_orn[e]'))=-1;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    if ~isempty(strmatch('SPMS',model.rxns,'exact'))
        model=changeObjective(model,'SPMS',1);
        %eFBA=optimizeCbModel(model,'min');
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'synthesis of spermidine from ornithine - SPMS';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    
    %% 295. synthesis of nitric oxide from arginine
    model=modelOri;
    model.lb(ismember(model.rxns,'EX_arg-L[e]'))=-1;model.ub(ismember(model.rxns,'EX_arg-L[e]'))=-1;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    if ~isempty(strmatch('NOS2',model.rxns,'exact'))
        model=changeObjective(model,'NOS2',1);
        %eFBA=optimizeCbModel(model,'min');
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'synthesis of nitric oxide from arginine - NOS2';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    
    %%  296. synthesis of cholesterol
    model=modelOri;
    for i = 1 : length(RPMI_composition)
        model = changeRxnBounds(model,RPMI_composition{i},-1,'l');
    end
    model=changeRxnBounds(model,'EX_o2[e]',-1,'l');
    if ~isempty(strmatch('DSREDUCr',model.rxns,'exact'))
        model=changeObjective(model,'DSREDUCr',1);
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'synthesis of cholesterol - DSREDUCr (with RPMI medium)';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    
    %% 297. denovo purine synthesis
    model=modelOri;
    model.lb(ismember(model.rxns,'EX_glc[e]'))=-1;model.ub(ismember(model.rxns,'EX_glc[e]'))=0;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    if ~isempty(strmatch('ADSL1',model.rxns,'exact'))
        model=changeObjective(model,'ADSL1',1);
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'de novo purine synthesis - ADSL1';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    
    %% 298. de novo purine synthesis - GMPS2
    model=modelOri;
    model.lb(ismember(model.rxns,'EX_glc[e]'))=-1;model.ub(ismember(model.rxns,'EX_glc[e]'))=0;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    if ~isempty(strmatch('GMPS2',model.rxns,'exact'))
        model=changeObjective(model,'GMPS2');
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'de novo purine synthesis - GMPS2';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    
    %% 299. salvage of purine bases
    model=modelOri;
    model.lb(ismember(model.rxns,'EX_glc[e]'))=-1;model.ub(ismember(model.rxns,'EX_glc[e]'))=0;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    if ~isempty(strmatch('ADPT',model.rxns,'exact'))
        model=changeObjective(model,'ADPT',1);
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'salvage of purine bases - ADPT';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    
    %% 300. salvage of purine bases - GUAPRT
    model=modelOri;
    model.lb(ismember(model.rxns,'EX_glc[e]'))=-1;model.ub(ismember(model.rxns,'EX_glc[e]'))=0;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    if ~isempty(strmatch('GUAPRT',model.rxns,'exact'))
        model=changeObjective(model,'GUAPRT',1);
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'salvage of purine bases - GUAPRT';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    
    %% 301. salvage of purine bases - HXPRT
    model=modelOri;
    model.lb(ismember(model.rxns,'EX_glc[e]'))=-1;model.ub(ismember(model.rxns,'EX_glc[e]'))=0;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    if ~isempty(strmatch('HXPRT',model.rxns,'exact'))
        model=changeObjective(model,'HXPRT',1);
        %eFBA=optimizeCbModel(model,'min');
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'salvage of purine bases - HXPRT';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    
    %% 302. purine catabolism
    model=modelOri;
    model.lb(ismember(model.rxns,'EX_glc[e]'))=-1;model.ub(ismember(model.rxns,'EX_glc[e]'))=0;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    if ~isempty(strmatch('XAOx',model.rxns,'exact'))
        model=changeObjective(model,'XAOx',1);
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'purine catabolism - XAOx';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    
    %% 303. pyrimidine synthesis (with hco3 uptake) - TMDS
    model=modelOri;
    model.lb(ismember(model.rxns,'EX_glc[e]'))=-1;model.ub(ismember(model.rxns,'EX_glc[e]'))=0;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    model=changeRxnBounds(model,'EX_hco3[e]',-1,'b');
    if ~isempty(strmatch('TMDS',model.rxns,'exact'))
        model=changeObjective(model,'TMDS',1);
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'pyrimidine synthesis (with hco3 uptake) - TMDS';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    
    %% 304. pyrimidine synthesis (with hco3 uptake) - CTPS2
    model=modelOri;
    model.lb(ismember(model.rxns,'EX_glc[e]'))=-1;model.ub(ismember(model.rxns,'EX_glc[e]'))=0;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    model=changeRxnBounds(model,'EX_hco3[e]',-1,'b');
    if ~isempty(strmatch('CTPS2',model.rxns,'exact'))
        model=changeObjective(model,'CTPS2',1);
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'pyrimidine synthesis (with hco3 uptake) - CTPS2';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    
    %% 305. pyrimidine catabolism
    model=modelOri;
    model.lb(ismember(model.rxns,'EX_glc[e]'))=-1;model.ub(ismember(model.rxns,'EX_glc[e]'))=0;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    model=changeRxnBounds(model,'EX_hco3[e]',-1,'b');
    if ~isempty(strmatch('UPPN',model.rxns,'exact'))
        model=changeObjective(model,'UPPN',1);
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'pyrimidine catabolism - UPPN';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    
    %% 306. 'pyrimidine catabolism - BUP2
    model=modelOri;
    model.lb(ismember(model.rxns,'EX_glc[e]'))=-1;model.ub(ismember(model.rxns,'EX_glc[e]'))=0;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    model=changeRxnBounds(model,'EX_hco3[e]',-1,'b');
    if ~isempty(strmatch('BUP2',model.rxns,'exact'))
        model=changeObjective(model,'BUP2',1);
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'pyrimidine catabolism - BUP2';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 307. fructose to glucose conversion
    model=modelOri;
    
    model.lb(ismember(model.rxns,'EX_fru[e]'))=-1;model.ub(ismember(model.rxns,'EX_fru[e]'))=-1;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    if ~isempty(strmatch('TRIOK',model.rxns,'exact'))
        model=changeObjective(model,'TRIOK',1);
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'fructose to glucose conversion - TRIOK';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 308. uptake and secretion of cholic acid
    model=modelOri;
    model=changeRxnBounds(model,'EX_o2[e]',-1,'l');
    model=changeRxnBounds(model,'EX_cholate[e]',-1,'l');
    model=changeRxnBounds(model,'EX_cholate[e]',1000,'u');
    % model=changeObjective(model,'CHOLATEt2u');
    %eFBA=optimizeCbModel(model,'min');
    %eFBA=optimizeCbModel(model);
    % TestSolution(k,1) = eFBA.f;
    % TestSolutionName{k,1} = 'uptake and secretion of cholic acid - CHOLATEt2u'; % SHOULD THIS BE MIN?
    % k = k +1;clear eFBA
    if ~isempty(strmatch('CHOLATEt3',model.rxns,'exact'))
        model=changeObjective(model,'CHOLATEt3',1);
        try
            model.osenseStr = 'min';
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;   
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'uptake of cholic acid - CHOLATEt3';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    
    %     if ~isempty(strmatch('CHOLATEt3',model.rxns,'exact'))
    %         eFBA=optimizeCbModel(model,'min');
    %         TestSolution(k,1) = eFBA.f;
    %     else
    %         TestSolution(k,1) = NaN;
    %     end
    %     TestSolutionName{k,1} = 'secretion of cholic acid - CHOLATEt3';
    %  if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    
    %% 309. Uptake and secretion of glycocholate
    model=modelOri;
    model=changeRxnBounds(model,'EX_o2[e]',-1,'l');
    model=changeRxnBounds(model,'EX_gchola[e]',-1,'l');
    model=changeRxnBounds(model,'EX_gchola[e]',1000,'u');
    % model=changeObjective(model,'GCHOLAt2u');
    %eFBA=optimizeCbModel(model,'min');
    %eFBA=optimizeCbModel(model);
    %TestSolution(k,1) = eFBA.f;
    %TestSolutionName{k,1} = 'uptake and secretion of cholic glycocholate - GCHOLAt2u';
    %k = k +1;clear eFBA
    if ~isempty(strmatch('GCHOLAt3',model.rxns,'exact'))
        model=changeObjective(model,'GCHOLAt3',1);
        %eFBA=optimizeCbModel(model,'min');
        try
            model.osenseStr = 'min';
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'uptake of cholic glycocholate - GCHOLAt3';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    
    %     if ~isempty(strmatch('GCHOLAt3',model.rxns,'exact'))
    %         eFBA=optimizeCbModel(model,'min');
    %         TestSolution(k,1) = eFBA.f;
    %     else
    %         TestSolution(k,1) = NaN;
    %     end
    %     TestSolutionName{k,1} = 'secretion of cholic glycocholate - GCHOLAt3';
    %  if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    
    %% 310. Uptake and secretion of tauro-cholate
    model=modelOri;
    model=changeRxnBounds(model,'EX_o2[e]',-1,'l');
    model=changeRxnBounds(model,'EX_tchola[e]',-1,'l');
    model=changeRxnBounds(model,'EX_tchola[e]',1000,'u');
    % model=changeObjective(model,'TCHOLAt2u');
    %eFBA=optimizeCbModel(model,'min');
    %eFBA=optimizeCbModel(model);
    %TestSolution(k,1) = eFBA.f;
    % TestSolutionName{k,1} = 'uptake and secretion of tauro-cholate - TCHOLAt2u';
    %k = k +1;clear eFBA
    if ~isempty(strmatch('TCHOLAt3',model.rxns,'exact'))
        model=changeObjective(model,'TCHOLAt3',1);
        % eFBA=optimizeCbModel(model,'min');
        try
            model.osenseStr = 'min';
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'uptake of tauro-cholate - TCHOLAt3';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %     if ~isempty(strmatch('TCHOLAt3',model.rxns,'exact'))
    %         eFBA=optimizeCbModel(model,'min');
    %         TestSolution(k,1) = eFBA.f;
    %     else
    %         TestSolution(k,1) = NaN;
    %     end
    %     TestSolutionName{k,1} = 'secretion of tauro-cholate - TCHOLAt3';
    %  if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    
    %% 311. Synthesis of fructose-6-phosphate from erythrose-4-phosphate (HMP shunt);
    model=modelOri;
    model.lb(ismember(model.rxns,'EX_glc[e]'))=-1;model.ub(ismember(model.rxns,'EX_glc[e]'))=0;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    if ~isempty(strmatch('TKT2',model.rxns,'exact'))
        model=changeObjective(model,'TKT2',1);
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Synthesis of fructose-6-phosphate from erythrose-4-phosphate (HMP shunt) - TKT2';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    
    %% 312. Malate to pyruvate (malic enzyme);
    model=modelOri;
    model.lb(ismember(model.rxns,'EX_glc[e]'))=-1;model.ub(ismember(model.rxns,'EX_glc[e]'))=0;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    if ~isempty(strmatch('ME2',model.rxns,'exact'))
        model=changeObjective(model,'ME2',1);
        eFBA = optimizeCbModel(model);
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Malate to pyruvate (malic enzyme) - ME2';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    
    %% 313. Malate to pyruvate (malic enzyme);
    model=modelOri;
    model.lb(ismember(model.rxns,'EX_glc[e]'))=-1;model.ub(ismember(model.rxns,'EX_glc[e]'))=0;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    if ~isempty(strmatch('ME2m',model.rxns,'exact'))
        model=changeObjective(model,'ME2m',1);
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Malate to pyruvate (malic enzyme) - ME2m';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    
    
    %% 314. Synthesis of urea (urea cycle);
    model=modelOri;
    for i = 1 : length(RPMI_composition)
        model = changeRxnBounds(model,RPMI_composition{i},-1,'l');
    end
    model=changeRxnBounds(model,'EX_o2[e]',-1,'l');
    if ~isempty(strmatch('ARGN',model.rxns,'exact'))
        model=changeObjective(model,'ARGN',1);
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Synthesis of urea (urea cycle) - ARGN (with RPMI medium)';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    
    %% 315. Cysteine to pyruvate
    model=modelOri;
    model=changeRxnBounds(model,'EX_o2[e]',-1,'l');
    model=changeRxnBounds(model,'EX_cys-L[e]',-1,'b');
    if ~isempty(strmatch('3SPYRSP',model.rxns,'exact'))
        model=changeObjective(model,'3SPYRSP',1);
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Cysteine to pyruvate - 3SPYRSP';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    
    
    %% 316. Methionine to cysteine  (check for dependancy over pe_hs);
    model=modelOri;
    model=changeRxnBounds(model,'EX_met_L[e]',-1,'b');
    model=changeRxnBounds(model,'EX_o2[e]',-1,'l');
    model=changeRxnBounds(model,'EX_pe_hs[e]',-1,'l');
    if ~isempty(strmatch('CYSTGL',model.rxns,'exact'))
        model=changeObjective(model,'CYSTGL',1);
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Methionine to cysteine - CYSTGL';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    
    %% 317. Synthesis of triacylglycerol (TAG reformation); (check for dependancy over dag_hs and RTOTAL3);
    model=modelOri;
    model=changeRxnBounds(model,'EX_o2[e]',-1,'l');
    model=changeRxnBounds(model,'EX_dag_hs[e]',-1,'l');
    model=changeRxnBounds(model,'EX_Rtotal3[e]',-1,'l');
    if ~isempty(strmatch('DGAT',model.rxns,'exact'))
        model=changeObjective(model,'DGAT');
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Synthesis of triacylglycerol (TAG reformation) - DGAT';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    
    %% 318. Phosphatidylcholine synthesis (check for dependancy over pe_hs);
    model=modelOri;
    model=changeRxnBounds(model,'EX_o2[e]',-1,'l');
    model=changeRxnBounds(model,'EX_pe_hs[e]',-1,'l');
    if ~isempty(strmatch('PETOHMm_hs',model.rxns,'exact'))
        model=changeObjective(model,'PETOHMm_hs',1);
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Phosphatidylcholine synthesis - PETOHMm_hs';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 319. Synthesis of FMN from riboflavin
    model=modelOri;
    model.lb(ismember(model.rxns,'EX_glc[e]'))=-1;model.ub(ismember(model.rxns,'EX_glc[e]'))=0;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=0;
    model=changeRxnBounds(model,'EX_ribflv[e]',-1,'b');
    if ~isempty(strmatch('RBFK',model.rxns,'exact'))
        model=changeObjective(model,'RBFK',1);
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Synthesis of FMN from riboflavin - RBFK';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    
    %% 320. synthesis of FAD from riboflavin
    model=modelOri;
    model.lb(ismember(model.rxns,'EX_glc[e]'))=-1;model.ub(ismember(model.rxns,'EX_glc[e]'))=0;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=0;
    model=changeRxnBounds(model,'EX_ribflv[e]',-1,'b');
    if ~isempty(strmatch('FMNAT',model.rxns,'exact'))
        model=changeObjective(model,'FMNAT',1);
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'synthesis of FAD from riboflavin - FMNAT';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 321. Synthesis of 5-methyl-tetrahydrofolate from folic acid
    model=modelOri;
    model=changeRxnBounds(model,'EX_o2[e]',-1,'l');
    model=changeRxnBounds(model,'EX_fol[e]',-1,'b');
    if ~isempty(strmatch('MTHFR3',model.rxns,'exact'))
        model=changeObjective(model,'MTHFR3',1);
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Synthesis of 5-methyl-tetrahydrofolate from folic acid - MTHFR3';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 322. Putriscine to GABA
    model=modelOri;
    model=changeRxnBounds(model,'EX_o2[e]',-1,'l');
    model=changeRxnBounds(model,'EX_ptrc[e]',-1,'b');
    if ~isempty(strmatch('ABUTD',model.rxns,'exact'))
        model=changeObjective(model,'ABUTD',1);
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Putriscine to GABA - ABUTD';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    
    %% 323. Superoxide dismutase
    model=modelOri;
    model.lb(ismember(model.rxns,'EX_glc[e]'))=-1;model.ub(ismember(model.rxns,'EX_glc[e]'))=0;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    if ~isempty(strmatch('SPODMm',model.rxns,'exact'))
        model=changeObjective(model,'SPODMm',1);
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Superoxide dismutase - SPODMm';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    
    %% 324. Availability of bicarbonate from Carbonic anhydrase reaction
    model=modelOri;
    model.lb(ismember(model.rxns,'EX_glc[e]'))=-1;model.ub(ismember(model.rxns,'EX_glc[e]'))=0;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    if ~isempty(strmatch('H2CO3Dm',model.rxns,'exact'))
        model=changeObjective(model,'H2CO3Dm',1);
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Availability of bicarbonate from Carbonic anhydrase reaction - H2CO3Dm';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    
    %% 325. Regeneration of citrate (TCA cycle);
    model=modelOri;
    model.lb(ismember(model.rxns,'EX_glc[e]'))=-1;model.ub(ismember(model.rxns,'EX_glc[e]'))=0;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    if ~isempty(strmatch('CSm',model.rxns,'exact'))
        model=changeObjective(model,'CSm',1);
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Regeneration of citrate (TCA cycle) - CSm';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    
    
    %% 326. Histidine to FIGLU
    model=modelOri;
    model.lb(find(ismember(model.rxns,'EX_his-L[e]')))=-1;
    model.ub(find(ismember(model.rxns,'EX_his-L[e]')))=-1;
    model=changeRxnBounds(model,'EX_o2[e]',-40,'l');
    model=changeRxnBounds(model,'EX_o2[e]',-1,'u');
    if ~isempty(strmatch('IZPN',model.rxns,'exact'))
        model=changeObjective(model,'IZPN',1);
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Histidine to FIGLU - IZPN';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    
    
    %% 327. binding of guar gum fiber to bile acids
    model=modelOri;
    model=changeRxnBounds(model,'EX_gum[e]',-1,'l');
    model=changeRxnBounds(model,'EX_gchola[e]',-1,'l');
    if ~isempty(strmatch('EX_gumgchol[e]',model.rxns,'exact'))
        model=changeObjective(model,'EX_gumgchol[e]',1);
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'binding of guar gum fiber to bile acids - EX_gumgchol[e]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    
    model=modelOri;
    model=changeRxnBounds(model,'EX_tchola[e]',-1,'l');
    model=changeRxnBounds(model,'EX_gum[e]',-1,'l');
    
    if ~isempty(strmatch('GUMTCHOLe',model.rxns,'exact'))
        model=changeObjective(model,'GUMTCHOLe',1);
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'binding of guar gum fiber to bile acids - GUMTCHOLe';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    
    model=modelOri;
    if ~isempty(strmatch('GUMDCHAe',model.rxns,'exact'))
        model=changeRxnBounds(model,'EX_dchac[e]',-1,'l');
        model=changeRxnBounds(model,'EX_gum[e]',-1,'l');
        model=changeObjective(model,'GUMDCHAe',1);
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'binding of guar gum fiber to bile acids - GUMDCHAe';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    
    %% 328. binding of psyllium fiber to bile acids
    model=modelOri;
    model=changeRxnBounds(model,'EX_psyl[e]',-1,'l');
    model=changeRxnBounds(model,'EX_gchola[e]',-1,'l');
    if ~isempty(strmatch('PSYGCHe',model.rxns,'exact'))
        model=changeObjective(model,'PSYGCHe',1);
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'binding of psyllium fiber to bile acids - PSYGCHe';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    
    
    model=modelOri;
    model=changeRxnBounds(model,'EX_psyl[e]',-1,'l');
    model=changeRxnBounds(model,'EX_tchola[e]',-1,'l');
    if ~isempty(strmatch('PSYTCHe',model.rxns,'exact'))
        model=changeObjective(model,'PSYTCHe',1);
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'binding of psyllium fiber to bile acids - PSYTCHe';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    
    model=modelOri;
    if ~isempty(strmatch('PSYTDECHe',model.rxns,'exact'))
        model=changeRxnBounds(model,'EX_tdechola[e]',-1,'l');
        model=changeRxnBounds(model,'EX_psyl[e]',-1,'l');
        model=changeObjective(model,'PSYTDECHe',1);
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'binding of psyllium fiber to bile acids - PSYTDECHe';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    
    %% 329. binding to beta glucan fibers to bile acids
    model=modelOri;
    model=changeRxnBounds(model,'EX_bglc[e]',-1,'l');
    model=changeRxnBounds(model,'EX_gchola[e]',-1,'l');
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    if ~isempty(strmatch('BGLUGCHe',model.rxns,'exact'))
        model=changeObjective(model,'BGLUGCHe',1);
        %eFBA=optimizeCbModel(model,'min');
        try
            model.osenseStr = 'min';
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'binding to beta glucan fibers to bile acids - BGLUGCHe';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    
    model=modelOri;
    model=changeRxnBounds(model,'EX_bglc[e]',-1,'l');
    model=changeRxnBounds(model,'EX_tchola[e]',-1,'l');
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    if ~isempty(strmatch('BGLUTCHLe',model.rxns,'exact'))
        model=changeObjective(model,'BGLUTCHLe',1);
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'binding to beta glucan fibers to bile acids - BGLUTCHLe';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    
    model=modelOri;
    model=changeRxnBounds(model,'EX_bglc[e]',-1,'l');
    model=changeRxnBounds(model,'EX_tdechola[e]',-1,'l');
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    if ~isempty(strmatch('BGLUTDECHOe',model.rxns,'exact'))
        model=changeObjective(model,'BGLUTDECHOe',1);
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'binding to beta glucan fibers to bile acids - BGLUTDECHOe';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    
    %% 330. binding of pectin fiber to bile acids
    model=modelOri;
    model=changeRxnBounds(model,'EX_pect[e]',-1,'l');
    model=changeRxnBounds(model,'EX_gchola[e]',-1,'l');
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    if ~isempty(strmatch('PECGCHLe',model.rxns,'exact'))
        model=changeObjective(model,'PECGCHLe',1);
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'binding of pectin fiber to bile acids - PECGCHLe';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    
    model=modelOri;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    model=changeRxnBounds(model,'EX_pect[e]',-1,'l');
    model=changeRxnBounds(model,'EX_tchola[e]',-1,'l');
    if ~isempty(strmatch('PECTCHLe',model.rxns,'exact'))
        model=changeObjective(model,'PECTCHLe',1);
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'binding of pectin fiber to bile acids - PECTCHLe';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    
    model=modelOri;
    if ~isempty(strmatch('PECDCHe',model.rxns,'exact'))
        model=changeRxnBounds(model,'EX_dchac[e]',-1,'l');
        model=changeRxnBounds(model,'EX_pect[e]',-1,'l');
        model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
        model=changeObjective(model,'PECDCHe',1);
        eFBA = optimizeCbModel(model);
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'binding of pectin fiber to bile acids - PECDCHe';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    
    %% 331. heme synthesis
    model=modelOri;
    model.lb(ismember(model.rxns,'EX_glc[e]'))=-1;model.ub(ismember(model.rxns,'EX_glc[e]'))=0;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    if ~isempty(strmatch('FCLTm',model.rxns,'exact'))
        model=changeObjective(model,'FCLTm',1);
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'heme synthesis - FCLTm';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    
    %% 332. heme degradation
    model=modelOri;
    model.lb(ismember(model.rxns,'EX_pheme[e]'))=-1;model.ub(ismember(model.rxns,'EX_pheme[e]'))=-1;
    model.lb(ismember(model.rxns,'EX_glc[e]'))=-1;model.ub(ismember(model.rxns,'EX_glc[e]'))=0;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=0;
    if ~isempty(strmatch('HOXG',model.rxns,'exact'))
        model=changeObjective(model,'HOXG',1);
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'heme degradation - HOXG';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    
end

%% these functions are new based on muscle and kidney work of SS
if strcmp(test,'all')|| strcmp(test,'Harvey')
    %% 333. Muscle objectives: valine -> pyruvate
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_val_L[e]'))=-1;model.ub(ismember(model.rxns,'EX_val_L[e]'))=-1;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    [model] = addSinkReactions(model,{'pyr[m]'},[0 100]);
    model.c(ismember(model.rxns,'sink_pyr[m]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'valine -> pyruvate';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 334. leucine -> pyruvate
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_leu_L[e]'))=-1;model.ub(ismember(model.rxns,'EX_leu_L[e]'))=-1;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    [model] = addSinkReactions(model,{'pyr[c]'},[0 100]);
    model.c(ismember(model.rxns,'sink_pyr[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'leucine -> pyruvate';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 335. isoleucine -> pyruvate
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_glc[e]'))=0;model.ub(ismember(model.rxns,'EX_glc[e]'))=0;
    model.lb(ismember(model.rxns,'EX_ile_L[e]'))=-1;model.ub(ismember(model.rxns,'EX_ile_L[e]'))=-1;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    [model] = addSinkReactions(model,{'pyr[c]'},[0 100]);
    model.c(ismember(model.rxns,'sink_pyr[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'isoleucine -> pyruvate';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 336. threonine -> alanine
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_glc[e]'))=0;model.ub(ismember(model.rxns,'EX_glc[e]'))=0;
    model.lb(ismember(model.rxns,'EX_thr_L[e]'))=-1;model.ub(ismember(model.rxns,'EX_thr_L[e]'))=-1;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=0;
    [model] = addSinkReactions(model,{'ala-L[c]'},[0 100]);
    model.c(ismember(model.rxns,'sink_ala-L[c]'))=1;
    model.c(ismember(model.rxns,'sink_ala_L[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'threonine -> alanine';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 337. aspartate -> pyruvate
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_glc[e]'))=0;model.ub(ismember(model.rxns,'EX_glc[e]'))=0;
    model.lb(ismember(model.rxns,'EX_asp_L[e]'))=-1;model.ub(ismember(model.rxns,'EX_asp_L[e]'))=-1;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    [model] = addSinkReactions(model,{'pyr[c]'},[0 100]);
    model.c(ismember(model.rxns,'sink_pyr[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'aspartate -> pyruvate';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 338. serine -> alanine
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_glc[e]'))=0;model.ub(ismember(model.rxns,'EX_glc[e]'))=0;
    model.lb(ismember(model.rxns,'EX_ser_L[e]'))=-1;model.ub(ismember(model.rxns,'EX_ser_L[e]'))=-1;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    [model] = addSinkReactions(model,{'ala-L[c]'},[0 100]);
    model.c(ismember(model.rxns,'sink_ala-L[c]'))=1;
    model.c(ismember(model.rxns,'sink_ala_L[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'serine -> alanine';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 339. glycine -> alanine
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_glc[e]'))=0;model.ub(ismember(model.rxns,'EX_glc[e]'))=0;
    model.lb(ismember(model.rxns,'EX_gly[e]'))=-1;model.ub(ismember(model.rxns,'EX_gly[e]'))=-1;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    [model] = addSinkReactions(model,{'ala-L[c]'},[0 100]);
    model.c(ismember(model.rxns,'sink_ala-L[c]'))=1;
    model.c(ismember(model.rxns,'sink_ala_L[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'glycine -> alanine';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 340. aspartate -> alanine
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_glc[e]'))=0;model.ub(ismember(model.rxns,'EX_glc[e]'))=0;
    model.lb(ismember(model.rxns,'EX_asp_L[e]'))=-1;model.ub(ismember(model.rxns,'EX_asp_L[e]'))=-1;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    [model] = addSinkReactions(model,{'ala-L[c]'},[0 100]);
    model.c(ismember(model.rxns,'sink_ala-L[c]'))=1;
    model.c(ismember(model.rxns,'sink_ala_L[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'aspartate -> alanine';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 341. tyrosine -> glutamine
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_glc[e]'))=0;model.ub(ismember(model.rxns,'EX_glc[e]'))=0;
    model.lb(ismember(model.rxns,'EX_tyr_L[e]'))=-1;model.ub(ismember(model.rxns,'EX_tyr_L[e]'))=-1;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    [model] = addSinkReactions(model,{'gln-L[c]'},[0 100]);
    model.c(ismember(model.rxns,'sink_gln-L[c]'))=1;
    model.c(ismember(model.rxns,'sink_gln_L[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'tyrosine -> glutamine';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 342. lysine -> glutamine
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_glc[e]'))=0;model.ub(ismember(model.rxns,'EX_glc[e]'))=0;
    model.lb(ismember(model.rxns,'EX_lys-L[e]'))=-1;model.ub(ismember(model.rxns,'EX_lys-L[e]'))=-1;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    [model] = addSinkReactions(model,{'gln-L[c]'},[0 100]);
    model.c(ismember(model.rxns,'sink_gln-L[c]'))=1;
    model.c(ismember(model.rxns,'sink_gln_L[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'lysine -> glutamine';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 343. phenylalanine -> glutamine
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_glc[e]'))=0;model.ub(ismember(model.rxns,'EX_glc[e]'))=0;
    model.lb(ismember(model.rxns,'EX_phe_L[e]'))=-1;model.ub(ismember(model.rxns,'EX_phe_L[e]'))=-1;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    [model] = addSinkReactions(model,{'gln-L[c]'},[0 100]);
    model.c(ismember(model.rxns,'sink_gln-L[c]'))=1;
    model.c(ismember(model.rxns,'sink_gln_L[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'phenylalanine -> glutamine';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 344. cysteine -> glutamine
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_glc[e]'))=0;model.ub(ismember(model.rxns,'EX_glc[e]'))=0;
    model.lb(ismember(model.rxns,'EX_cys-L[e]'))=-1;model.ub(ismember(model.rxns,'EX_cys-L[e]'))=-1;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    [model] = addSinkReactions(model,{'gln-L[c]'},[0 100]);
    model.c(ismember(model.rxns,'sink_gln-L[c]'))=1;
    model.c(ismember(model.rxns,'sink_gln_L[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'cysteine -> glutamine';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 345. cysteine -> alanine
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_glc[e]'))=0;model.ub(ismember(model.rxns,'EX_glc[e]'))=0;
    model.lb(ismember(model.rxns,'EX_cys-L[e]'))=-1;model.ub(ismember(model.rxns,'EX_cys-L[e]'))=-1;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    [model] = addSinkReactions(model,{'ala-L[c]'},[0 100]);
    model.c(ismember(model.rxns,'sink_ala-L[c]'))=1;
    model.c(ismember(model.rxns,'sink_ala_L[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'cysteine -> alanine';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 346. leucine -> glutamine
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_glc[e]'))=0;model.ub(ismember(model.rxns,'EX_glc[e]'))=0;
    model.lb(ismember(model.rxns,'EX_leu_L[e]'))=-1;model.ub(ismember(model.rxns,'EX_leu_L[e]'))=-1;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    [model] = addSinkReactions(model,{'gln-L[c]'},[0 100]);
    model.c(ismember(model.rxns,'sink_gln-L[c]'))=1;
    model.c(ismember(model.rxns,'sink_gln_L[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'leucine -> glutamine';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 347. leucine -> alanine
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_glc[e]'))=0;model.ub(ismember(model.rxns,'EX_glc[e]'))=0;
    model.lb(ismember(model.rxns,'EX_leu_L[e]'))=-1;model.ub(ismember(model.rxns,'EX_leu_L[e]'))=-1;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    [model] = addSinkReactions(model,{'ala-L[c]'},[0 100]);
    model.c(ismember(model.rxns,'sink_ala-L[c]'))=1;
    model.c(ismember(model.rxns,'sink_ala_L[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'leucine -> alanine';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 348. valine -> glutamine
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_glc[e]'))=0;model.ub(ismember(model.rxns,'EX_glc[e]'))=0;
    model.lb(ismember(model.rxns,'EX_val_L[e]'))=-1;model.ub(ismember(model.rxns,'EX_val_L[e]'))=-1;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    [model] = addSinkReactions(model,{'gln-L[c]'},[0 100]);
    model.c(ismember(model.rxns,'sink_gln-L[c]'))=1;
    model.c(ismember(model.rxns,'sink_gln_L[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'valine -> glutamine';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 349. valine -> alanine
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_glc[e]'))=0;model.ub(ismember(model.rxns,'EX_glc[e]'))=0;
    model.lb(ismember(model.rxns,'EX_val_L[e]'))=-1;model.ub(ismember(model.rxns,'EX_val_L[e]'))=-1;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    [model] = addSinkReactions(model,{'ala-L[c]'},[0 100]);
    model.c(ismember(model.rxns,'sink_ala-L[c]'))=1;
    model.c(ismember(model.rxns,'sink_ala_L[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'valine -> alanine';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 350. isoleucine -> glutamine
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_glc[e]'))=0;model.ub(ismember(model.rxns,'EX_glc[e]'))=0;
    model.lb(ismember(model.rxns,'EX_ile_L[e]'))=-1;model.ub(ismember(model.rxns,'EX_ile_L[e]'))=-1;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    [model] = addSinkReactions(model,{'gln-L[c]'},[0 100]);
    model.c(ismember(model.rxns,'sink_gln-L[c]'))=1;
    model.c(ismember(model.rxns,'sink_gln_L[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'isoleucine -> glutamine';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 351. isoleucine -> alanine
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_glc[e]'))=0;model.ub(ismember(model.rxns,'EX_glc[e]'))=0;
    model.lb(ismember(model.rxns,'EX_ile_L[e]'))=-1;model.ub(ismember(model.rxns,'EX_ile_L[e]'))=-1;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    [model] = addSinkReactions(model,{'ala-L[c]'},[0 100]);
    model.c(ismember(model.rxns,'sink_ala-L[c]'))=1;
    model.c(ismember(model.rxns,'sink_ala_L[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'isoleucine -> alanine';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 352. methionine -> glutamine
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_glc[e]'))=0;model.ub(ismember(model.rxns,'EX_glc[e]'))=0;
    model.lb(ismember(model.rxns,'EX_met_L[e]'))=-1;model.ub(ismember(model.rxns,'EX_met_L[e]'))=-1;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    [model] = addSinkReactions(model,{'gln-L[c]'},[0 100]);
    model.c(ismember(model.rxns,'sink_gln-L[c]'))=1;
    model.c(ismember(model.rxns,'sink_gln_L[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'methionine -> glutamine';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 353. methionine -> alanine
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_glc[e]'))=0;model.ub(ismember(model.rxns,'EX_glc[e]'))=0;
    model.lb(ismember(model.rxns,'EX_met_L[e]'))=-1;model.ub(ismember(model.rxns,'EX_met_L[e]'))=-1;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    [model] = addSinkReactions(model,{'ala-L[c]'},[0 100]);
    model.c(ismember(model.rxns,'sink_ala-L[c]'))=1;
    model.c(ismember(model.rxns,'sink_ala_L[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'methionine -> alanine';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 354. arginine -> ornithine
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_glc[e]'))=0;model.ub(ismember(model.rxns,'EX_glc[e]'))=0;
    model.lb(ismember(model.rxns,'EX_arg-L[e]'))=-1;model.ub(ismember(model.rxns,'EX_arg-L[e]'))=-1;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    [model] = addSinkReactions(model,{'orn[c]'},[0 100]);
    model.c(ismember(model.rxns,'sink_orn[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'arginine -> ornithine';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 355. arginine -> proline
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_glc[e]'))=0;model.ub(ismember(model.rxns,'EX_glc[e]'))=0;
    model.lb(ismember(model.rxns,'EX_arg-L[e]'))=-1;model.ub(ismember(model.rxns,'EX_arg-L[e]'))=-1;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    [model] = addSinkReactions(model,{'pro-L[c]'},[0 100]);
    model.c(ismember(model.rxns,'sink_pro-L[c]'))=1;
    model.c(ismember(model.rxns,'sink_pro_L[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'arginine -> proline';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 356. ornithine -> putrescine
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_glc[e]'))=0;model.ub(ismember(model.rxns,'EX_glc[e]'))=0;
    model.lb(ismember(model.rxns,'EX_orn[e]'))=-1;model.ub(ismember(model.rxns,'EX_orn[e]'))=-1;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    [model] = addSinkReactions(model,{'ptrc[c]'},[0 100]);
    model.c(ismember(model.rxns,'sink_ptrc[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'ornithine -> putrescine';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 357. glutamate -> glutamine
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_glc[e]'))=0;model.ub(ismember(model.rxns,'EX_glc[e]'))=0;
    model.lb(ismember(model.rxns,'EX_glu-L[e]'))=-1;model.ub(ismember(model.rxns,'EX_glu-L[e]'))=-1;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    [model] = addSinkReactions(model,{'gln-L[c]'},[0 100]);
    model.c(ismember(model.rxns,'sink_gln-L[c]'))=1;
    model.c(ismember(model.rxns,'sink_gln_L[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'glutamate -> glutamine';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 358. methionine -> spermine
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_glc[e]'))=0;model.ub(ismember(model.rxns,'EX_glc[e]'))=0;
    model.lb(ismember(model.rxns,'EX_met_L[e]'))=-1;model.ub(ismember(model.rxns,'EX_met_L[e]'))=-1;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    [model] = addSinkReactions(model,{'sprm[c]'},[0 100]);
    model.c(ismember(model.rxns,'sink_sprm[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'methionine -> spermine';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 359. methionine -> spermidine
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_glc[e]'))=0;model.ub(ismember(model.rxns,'EX_glc[e]'))=0;
    model.lb(ismember(model.rxns,'EX_met_L[e]'))=-1;model.ub(ismember(model.rxns,'EX_met_L[e]'))=-1;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    [model] = addSinkReactions(model,{'spmd[c]'},[0 100]);
    model.c(ismember(model.rxns,'sink_spmd[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'methionine -> spermidine';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 360. spermidine -> putrescine
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_glc[e]'))=0;model.ub(ismember(model.rxns,'EX_glc[e]'))=0;
    model.lb(ismember(model.rxns,'EX_spmd[e]'))=-1;model.ub(ismember(model.rxns,'EX_spmd[e]'))=-1;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    [model] = addSinkReactions(model,{'ptrc[c]'},[0 100]);
    model.c(ismember(model.rxns,'sink_ptrc[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'spermidine -> putrescine';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 361. ADP -> ATP/ adenylate kinase
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    model.lb(ismember(model.rxns,'EX_glc[e]'))=-1;model.ub(ismember(model.rxns,'EX_glc[e]'))=0;
    if ~isempty(strmatch('AK1',model.rxns,'exact'))
        model.c(ismember(model.rxns,'AK1'))=1;
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'ADP -> ATP/ adenylate kinase';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 362. ADP -> ATP/ adenylate kinase
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    model.lb(ismember(model.rxns,'EX_glc[e]'))=-1;model.ub(ismember(model.rxns,'EX_glc[e]'))=0;
    if ~isempty(strmatch('AK1',model.rxns,'exact'))
        model.c(ismember(model.rxns,'AK1m'))=1;
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'ADP -> ATP/ adenylate kinase (mitochondrial)';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 363. phosphocreatine -> creatine/ cytosolic creatine kinase
    model = modelOri;
    model.c(find(model.c)) = 0;
    model = addReaction(model,'EX_pcreat[e]','pcreat[e] <=>');
    model.lb(ismember(model.rxns,'EX_glc[e]'))=-1;model.ub(ismember(model.rxns,'EX_glc[e]'))=0;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=0;
    model.lb(ismember(model.rxns,'EX_pcreat[e]'))=-1;model.ub(ismember(model.rxns,'EX_pcreat[e]'))=-1;
    [model] = addSinkReactions(model,{'creat[c]'},[0 100]);
    model.c(ismember(model.rxns,'sink_creat[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'phosphocreatine -> creatine/ cytosolic creatine kinase';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 364. creatine -> phosphocreatine/mitochondrial creatine kinase
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_creat[e]'))=-1;model.ub(ismember(model.rxns,'EX_creat[e]'))=-1;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    [model] = addSinkReactions(model,{'pcreat[c]'},[0 100]);
    model.c(ismember(model.rxns,'sink_pcreat[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'creatine -> phosphocreatine/mitochondrial creatine kinase';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 365. fructose -> lactate/ oxidation of fructose
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_fru[e]'))=-1;model.ub(ismember(model.rxns,'EX_fru[e]'))=-1;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=0;
    [model] = addSinkReactions(model,{'lac-L[c]'},[0 100]);
    model.c(ismember(model.rxns,'sink_lac-L[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'fructose -> lactate/ oxidation of fructose';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 366. fructose -> glycogen/ glycogenesis
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_fru[e]'))=-1;model.ub(ismember(model.rxns,'EX_fru[e]'))=-1;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=0;
    [model] = addSinkReactions(model,{'glygn2[c]'},[0 100]);
    model.c(ismember(model.rxns,'sink_glygn2[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'fructose -> glycogen/ glycogenesis';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 367. glucose -> erythrose/ HMP shunt
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_glc[e]'))=-1;model.ub(ismember(model.rxns,'EX_glc[e]'))=-1;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    [model] = addSinkReactions(model,{'e4p[c]'},[0 100]);
    model.c(ismember(model.rxns,'sink_e4p[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'glucose -> erythrose/ HMP shunt';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 368. tag_hs[c] -> mag_hs[c]/ lipolysis
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_glc[e]'))=0;model.ub(ismember(model.rxns,'EX_glc[e]'))=0;
    model.lb(ismember(model.rxns,'EX_tag_hs[e]'))=-1;model.ub(ismember(model.rxns,'EX_tag_hs[e]'))=-1;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=0;
    [model] = addSinkReactions(model,{'mag-hs[c]'},[0 100]);
    model.c(ismember(model.rxns,'sink_mag-hs[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'tag_hs[c] -> mag_hs[c]/ lipolysis';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 369. tag_hs[c] -> glyc[c]/ lipolysis
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_glc[e]'))=0;model.ub(ismember(model.rxns,'EX_glc[e]'))=0;
    model.lb(ismember(model.rxns,'EX_tag_hs[e]'))=-1;model.ub(ismember(model.rxns,'EX_tag_hs[e]'))=-1;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    [model] = addSinkReactions(model,{'glyc[c]'},[0 100]);
    model.c(ismember(model.rxns,'sink_glyc[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'tag_hs[c] -> glyc[c]/ lipolysis';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 370. pmtcoa -> acetylCoA/ beta oxidation from pmtcoa
    model = modelOri;
    %         for i = 1 : length(RPMI_composition)
    %         model = changeRxnBounds(model,RPMI_composition{i},-1,'l');
    %     end
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_hdca[e]'))=-1;model.ub(ismember(model.rxns,'EX_hdca[e]'))=-1;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    model.lb(find(ismember(model.rxns,'sink_coa[c]')))=-1;
    model.ub(find(ismember(model.rxns,'sink_coa[c]')))=1;
    [model] = addSinkReactions(model,{'accoa[m]'},[0 100]);
    model.c(ismember(model.rxns,'sink_accoa[m]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'pmtcoa -> acetylCoA/ beta oxidation from pmtcoa';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 371. odecoa -> acetylCoA/ beta oxidation from oleic acid
    model = modelOri;
    %         for i = 1 : length(RPMI_composition)
    %         model = changeRxnBounds(model,RPMI_composition{i},-1,'l');
    %     end
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_ocdcea[e]'))=-1;model.ub(ismember(model.rxns,'EX_ocdcea[e]'))=-1;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=0;
    model.lb(find(ismember(model.rxns,'sink_coa[c]')))=-1;
    model.ub(find(ismember(model.rxns,'sink_coa[c]')))=1;
    [model] = addSinkReactions(model,{'accoa[m]'},[0 100]);
    model.c(ismember(model.rxns,'sink_accoa[m]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'odecoa -> acetylCoA/ beta oxidation from oleic acid (with RPMI medium)';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 372. lnlccoa -> acetylCoA/ beta oxidation from linoleic acid
    model = modelOri;
    %         for i = 1 : length(RPMI_composition)
    %         model = changeRxnBounds(model,RPMI_composition{i},-1,'l');
    %     end
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_lnlc[e]'))=-1;model.ub(ismember(model.rxns,'EX_lnlc[e]'))=-1;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=0;
    model.lb(find(ismember(model.rxns,'sink_coa[c]')))=-1;
    model.ub(find(ismember(model.rxns,'sink_coa[c]')))=1;
    [model] = addSinkReactions(model,{'accoa[m]'},[0 100]);
    model.c(ismember(model.rxns,'sink_accoa[m]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'lnlccoa -> acetylCoA/ beta oxidation from linoleic acid (with RPMI medium)';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 373. glycerol -> dhap/ glycerol utilizing machinery
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_glc[e]'))=0;model.ub(ismember(model.rxns,'EX_glc[e]'))=0;
    model.lb(ismember(model.rxns,'EX_glyc[e]'))=-1;model.ub(ismember(model.rxns,'EX_glyc[e]'))=-1;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    [model] = addSinkReactions(model,{'dhap[c]'},[0 100]);
    model.c(ismember(model.rxns,'sink_dhap[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'glycerol -> dhap/ glycerol utilizing machinery';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 374. adenine -> amp/ salvage of adenine
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_glc[e]'))=0;model.ub(ismember(model.rxns,'EX_glc[e]'))=0;
    model.lb(ismember(model.rxns,'EX_adn[e]'))=-1;model.ub(ismember(model.rxns,'EX_adn[e]'))=-1;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    [model] = addSinkReactions(model,{'amp[c]'},[0 100]);
    model.c(ismember(model.rxns,'sink_amp[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'adenine -> amp/ salvage of adenine';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 375. hypoxanthine -> imp/ salvage of hypoxanthine
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_hxan[e]'))=-1;model.ub(ismember(model.rxns,'EX_hxan[e]'))=-1;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    [model] = addSinkReactions(model,{'imp[c]'},[0 100]);
    model.c(ismember(model.rxns,'INSK'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'hypoxanthine -> imp/ salvage of hypoxanthine';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 376. guanine -> gmp/ salvage of guanine
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_gua[e]'))=-1;model.ub(ismember(model.rxns,'EX_gua[e]'))=-1;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=0;
    [model] = addSinkReactions(model,{'prpp[c]','gmp[c]'},[-1 0;0 100]);
    model.c(ismember(model.rxns,'GUAPRT'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'guanine -> gmp/ salvage of guanine';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 377. ribose -> imp/ denovo purine synthesis
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_glc[e]'))=0;model.ub(ismember(model.rxns,'EX_glc[e]'))=0;
    model.lb(ismember(model.rxns,'EX_rib_D[e]'))=-1;model.ub(ismember(model.rxns,'EX_rib_D[e]'))=-1;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    [model] = addSinkReactions(model,{'imp[c]'},[0 100]);
    model.c(ismember(model.rxns,'sink_imp[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'ribose -> imp/ denovo purine synthesis';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 378. thymd -> thym/ thymidine phosphorylase
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_glc[e]'))=0;model.ub(ismember(model.rxns,'EX_glc[e]'))=0;
    model.lb(ismember(model.rxns,'EX_thymd[e]'))=-1;model.ub(ismember(model.rxns,'EX_thymd[e]'))=-1;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    [model] = addSinkReactions(model,{'thym[c]'},[0 100]);
    model.c(ismember(model.rxns,'sink_thym[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'thymd -> thym/ thymidine phosphorylase';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 379. glutamine -> cmp/ pyrimidine synthesis
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_glc[e]'))=0;model.ub(ismember(model.rxns,'EX_glc[e]'))=0;
    model.lb(ismember(model.rxns,'EX_gln-L[e]'))=-1;model.ub(ismember(model.rxns,'EX_gln-L[e]'))=-1;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    [model] = addSinkReactions(model,{'cmp[c]'},[0 100]);
    model.c(ismember(model.rxns,'sink_cmp[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'glutamine -> cmp/ pyrimidine synthesis';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 380. glutamine -> dtmp/ pyrimidine synthesis
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_glc[e]'))=0;model.ub(ismember(model.rxns,'EX_glc[e]'))=0;
    model.lb(ismember(model.rxns,'EX_gln-L[e]'))=-1;model.ub(ismember(model.rxns,'EX_gln-L[e]'))=-1;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    [model] = addSinkReactions(model,{'dtmp[c]'},[0 100]);
    model.c(ismember(model.rxns,'sink_dtmp[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'glutamine -> dtmp/ pyrimidine synthesis';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 381. Kidney objectives: citr_L[c] -> arg_L[c]
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'citr-L[c]','arg-L[c]'},[-1 -1; 0 100]);
    model.lb(find(ismember(model.rxns,'sink_citr[c]')))=-1;
    model.ub(find(ismember(model.rxns,'sink_citr[c]')))=-1;
    model.c(ismember(model.rxns,'sink_arg-L[c]'))=1;
    model.c(ismember(model.rxns,'sink_arg_L[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'citr_L[c] -> arg_L[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 382. cys_L[c] -> taur[c]
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'cys-L[c]','taur[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_taur[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'cys_L[c] -> taur[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 383. gly[c] -> orn[c]
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'gly[c]','orn[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_orn[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'gly[c] -> orn[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 384. citr_L[c] -> urea[c]/ partial urea cycle in kidney
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'citr-L[c]','urea[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_urea[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'citr_L[c] -> urea[c]/ partial urea cycle in kidney';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 385. gthrd[c] -> glycine[c]/ glutathione breakdown via ?-glutamyl-transeptidase
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_glc[e]'))=-1;model.ub(ismember(model.rxns,'EX_glc[e]'))=0;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=0;
    [model] = addSinkReactions(model,{'gly[c]','gthrd[c]'},[-1 -1; 0 100]);
    model.lb(find(ismember(model.rxns,'sink_gly[c]')))=-1;
    model.ub(find(ismember(model.rxns,'sink_gly[c]')))=-1;
    model.c(ismember(model.rxns,'sink_gthrd[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'gthrd[c] -> glycine[c]/ glutathione breakdown via glutamyl-transeptidase';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 386. pro_L[c] -> GABA[c]/ GABA synthesis in kidney
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'pro-L[c]','4abut[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_4abut[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'pro_L[c] -> GABA[c]/ GABA synthesis in kidney';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 387. pro_L[c] -> orn[c]
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'pro-L[c]','orn[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_orn[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'pro_L[c] -> orn[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 388. met_L[c] -> hcys_L[c]
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'met-L[c]','hcys-L[c]'},[-1 -1; 0 100]);
    model.lb(find(ismember(model.rxns,'sink_met_L[c]')))=-1;
    model.ub(find(ismember(model.rxns,'sink_met_L[c]')))=-1;
    model.c(ismember(model.rxns,'sink_hcys-L[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'met_L[c] -> hcys_L[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 389. hcys_L[c] -> met_L[c]
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'hcys-L[c]','met-L[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_met-L[c]'))=1;
    model.c(ismember(model.rxns,'sink_met_L[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'hcys_L[c] -> met_L[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 390. hcys_L[c] -> cys_L[c]
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=0;
    [model] = addSinkReactions(model,{'hcys-L[c]','cys-L[c]'},[-1 -1; 0 100]);
    model.lb(find(ismember(model.rxns,'sink_ser_L[c]')))=-1;
    model.ub(find(ismember(model.rxns,'sink_ser_L[c]')))=-1;
    model.c(ismember(model.rxns,'sink_cys-L[c]'))=1;
    model.c(ismember(model.rxns,'sink_cys_L[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'hcys_L[c] -> cys_L[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 391. 'lys-L[c] -> glu_L[c] / lysine degradation
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'lys-L[c]','glu-L[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_glu-L[c]'))=1;
    model.c(ismember(model.rxns,'sink_glu_L[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'lys-L[c] -> glu_L[c] / lysine degradation';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 392. trp-L[c] -> trypta[c] / tryptophan degradation
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'trp-L[c]','trypta[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_trypta[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'trp-L[c] -> trypta[c] / tryptophan degradation';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 393. kynate[c] -> nicotinamide[c] / nicotinamide from tryptophan metabolite
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'kynate[c]','nicrnt[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_nicrnt[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'kynate[c] -> nicotinamide[c] / nicotinamide from tryptophan metabolite';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 394. pyr[c] -> lac-L[c]/ lactate dehydrogenase
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'pyr[c]','lac-L[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_lac-L[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'pyr[c] -> lac-L[c]/ lactate dehydrogenase';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 395. ATP max, aerobic, pyruvate/ pyruvate dehydrogenase-->TCA->energy
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_pyr[e]'))=-1;model.ub(ismember(model.rxns,'EX_pyr[e]'))=-1;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    
    if ~isempty(strmatch('ATPM',model.rxns,'exact'))
        model.c(ismember(model.rxns,'ATPM'))=1;
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'ATP max, aerobic, pyruvate/ pyruvate dehydrogenase-->TCA->energy';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 396. gal[c] -> udpg[c]/ galactose utilization
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'gal[c]','udpg[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_udpg[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'gal[c] -> udpg[c]/ galactose utilization';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 397. fru[c] -> lac_L[c]/ fructose conversion to glucose & utilization
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'fru[c]','lac-L[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_lac-L[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'fru[c] -> lac_L[c]/ fructose conversion to glucose & utilization';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 398. malcoa[c] -> eicostetcoa[c]/ fatty acid elongation
    model = modelOri;
    for i = 1 : length(RPMI_composition)
        model = changeRxnBounds(model,RPMI_composition{i},-1,'l');
    end
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'malcoa[c]','eicostetcoa[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_eicostetcoa[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'malcoa[c] -> eicostetcoa[c]/ fatty acid elongation (wtih RPMI medium)';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 399. accoa[c] -> chsterol[r]
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'accoa[c]','chsterol[r]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_chsterol[r]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'accoa[c] -> chsterol[r]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 400. inost[c] -> glac[r]
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'inost[c]','glac[r]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_glac[r]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'inost[c] -> glac[r]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 401. pail_hs[c] -> pail4p_hs[c]/ inositol kinase
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'pail_hs[c]','pail4p_hs[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_pail4p_hs[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'pail_hs[c] -> pail4p_hs[c]/ inositol kinase';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 402. arachd[c] -> prostgh2[c]/ prostaglandin synthesis
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'arachd[c]','prostgh2[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_prostgh2[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'arachd[c] -> prostgh2[c]/ prostaglandin synthesis';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 403. arachd[c] -> prostgd2[r]/ prostaglandin synthesis
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'arachd[c]','prostgd2[r]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_prostgd2[r]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'arachd[c] -> prostgd2[r]/ prostaglandin synthesis';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 404. arachd[c] -> prostge2[r]/ prostaglandin synthesis
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'arachd[c]','prostge2[r]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_prostge2[r]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'arachd[c] -> prostge2[r]/ prostaglandin synthesis';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 405. arachd[c] -> prostgi2[r]/ prostaglandin synthesis
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'arachd[c]','prostgi2[r]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_prostgi2[r]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'arachd[c] -> prostgi2[r]/ prostaglandin synthesis';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 406. 25hvitd3[m] -> 2425dhvitd3[m]/ 24,25-dihydroxycalciol synthesis
    model = modelOri;
    model.c(find(model.c)) = 0;
    for i = 1 : length(RPMI_composition)
        model = changeRxnBounds(model,RPMI_composition{i},-1,'l');
    end
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=0;
    [model] = addSinkReactions(model,{'25hvitd3[m]','2425dhvitd3[m]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_2425dhvitd3[m]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = '25hvitd3[m] -> 2425dhvitd3[m]/ 24,25-dihydroxycalciol synthesis (with RPMI medium)';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 407. caro[c] -> retinal[c]/ vitamin A synthesis
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'caro[c]','retinal[c]'},[-1 -1; 0 100]);
    model.c(ismember(model.rxns,'sink_retinal[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'caro[c] -> retinal[c]/ vitamin A synthesis';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
 
    %% 408. synthesis of glutamate from ornithine
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_orn[e]'))=-1;model.ub(ismember(model.rxns,'EX_orn[e]'))=-1;
    [model] = addDemandReaction(model,'glu-L[c]');
    model.c(ismember(model.rxns,'sink_retinal[c]'))=1;
    if find(model.c)>0
        try
            eFBA=entropicFluxBalanceAnalysis(model,param);
            eFBA.f=sum(eFBA.v(find(model.c>0)));
        catch
            eFBA.f=NaN;
        end
        TestSolution(k,1) = eFBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'synthesis of glutamate from ornithine';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 409. synthesis of proline from ornithine
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_orn[e]'))=-1;model.ub(ismember(model.rxns,'EX_orn[e]'))=-1;
    [model] = addDemandReaction(model,'pro-L[m]');
    model.c(ismember(model.rxns,'DM_pro-L[m]'))=1;
    try
        eFBA=entropicFluxBalanceAnalysis(model,param);
        eFBA.f=sum(eFBA.v(find(model.c>0)));
    catch
        eFBA.f=NaN;
    end
    TestSolution(k,1) = eFBA.f;
    TestSolutionName{k,1} = 'synthesis of proline from ornithine';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 410. visual cycle in retina
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'retinol-cis-11[c]','retinal[c]'},[-1 -1; 0 1000]);
    model.c(ismember(model.rxns,'sink_retinal[c]'))=1;
    eFBA = optimizeCbModel(model);
    TestSolution(k,1) = eFBA.f;
    TestSolutionName{k,1} = 'visual cycle in retina';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 411. pail_hs[c] -> pchol_hs[c]
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'pail_hs[c]','pchol-hs[c]'},[-1 -1; 0 1000]);
    model.c(ismember(model.rxns,'sink_pchol-hs[c]'))=1;
    try
        eFBA=entropicFluxBalanceAnalysis(model,param);
        eFBA.f=sum(eFBA.v(find(model.c>0)));
    catch
        eFBA.f=NaN;
    end
    TestSolution(k,1) = eFBA.f;
    TestSolutionName{k,1} = 'pail_hs[c] -> pchol_hs[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 412. pail_hs[c] -> pe_hs[c]
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'pail_hs[c]','pe_hs[c]'},[-1 -1; 0 1000]);
    model.c(ismember(model.rxns,'sink_pe_hs[c]'))=1;
    try
        eFBA=entropicFluxBalanceAnalysis(model,param);
        eFBA.f=sum(eFBA.v(find(model.c>0)));
    catch
        eFBA.f=NaN;
    end
    TestSolution(k,1) = eFBA.f;
    TestSolutionName{k,1} = 'pail_hs[c] -> pe_hs[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 413. pail_hs[c] -> ps_hs[c]
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'pail_hs[c]','ps-hs[c]'},[-1 -1; 0 1000]);
    model.c(ismember(model.rxns,'sink_ps-hs[c]'))=1;
    try
        eFBA=entropicFluxBalanceAnalysis(model,param);
        eFBA.f=sum(eFBA.v(find(model.c>0)));
    catch
        eFBA.f=NaN;
    end
    TestSolution(k,1) = eFBA.f;
    TestSolutionName{k,1} = 'pail_hs[c] -> ps_hs[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 414. pail_hs[c] -> g3pc[c]
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'pail_hs[c]','g3pc[c]'},[-1 -1; 0 1000]);
    model.c(ismember(model.rxns,'sink_g3pc[c]'))=1;
    try
        eFBA=entropicFluxBalanceAnalysis(model,param);
        eFBA.f=sum(eFBA.v(find(model.c>0)));
    catch
        eFBA.f=NaN;
    end
    TestSolution(k,1) = eFBA.f;
    TestSolutionName{k,1} = 'pail_hs[c] -> g3pc[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 415. dag_hs[c] -> pchol_hs[c]
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'dag_hs[c]','pchol-hs[c]'},[-1 -1; 0 1000]);
    model.c(ismember(model.rxns,'sink_pchol-hs[c]'))=1;
    try
        eFBA=entropicFluxBalanceAnalysis(model,param);
        eFBA.f=sum(eFBA.v(find(model.c>0)));
    catch
        eFBA.f=NaN;
    end
    TestSolution(k,1) = eFBA.f;
    TestSolutionName{k,1} = 'dag_hs[c] -> pchol_hs[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 416. dag_hs[c] -> pe_hs[c]
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'dag_hs[c]','pe_hs[c]'},[-1 -1; 0 1000]);
    model.c(ismember(model.rxns,'sink_pe_hs[c]'))=1;
    try
        eFBA=entropicFluxBalanceAnalysis(model,param);
        eFBA.f=sum(eFBA.v(find(model.c>0)));
    catch
        eFBA.f=NaN;
    end
    TestSolution(k,1) = eFBA.f;
    TestSolutionName{k,1} = 'dag_hs[c] -> pe_hs[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 417. dag_hs[c] -> clpn_hs[c]
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'dag_hs[c]','clpn-hs[c]'},[-1 -1; 0 1000]);
    model.c(ismember(model.rxns,'sink_clpn-hs[c]'))=1;
    try
        eFBA=entropicFluxBalanceAnalysis(model,param);
        eFBA.f=sum(eFBA.v(find(model.c>0)));
    catch
        eFBA.f=NaN;
    end
    TestSolution(k,1) = eFBA.f;
    TestSolutionName{k,1} = 'dag_hs[c] -> clpn_hs[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 418. dag_hs[c] -> pgp_hs[c]
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'dag_hs[c]','pgp-hs[c]'},[-1 -1; 0 1000]);
    model.c(ismember(model.rxns,'sink_pgp-hs[c]'))=1;
    try
        eFBA=entropicFluxBalanceAnalysis(model,param);
        eFBA.f=sum(eFBA.v(find(model.c>0)));
    catch
        eFBA.f=NaN;
    end
    TestSolution(k,1) = eFBA.f;
    TestSolutionName{k,1} = 'dag_hs[c] -> pgp_hs[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 419. bhb[m] -> acac[m]/ ketone body utilization
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'bhb[m]','acac[m]'},[-1 -1; 0 1000]);
    model.c(ismember(model.rxns,'sink_acac[m]'))=1;
    try
        eFBA=entropicFluxBalanceAnalysis(model,param);
        eFBA.f=sum(eFBA.v(find(model.c>0)));
    catch
        eFBA.f=NaN;
    end
    TestSolution(k,1) = eFBA.f;
    TestSolutionName{k,1} = 'bhb[m] -> acac[m]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 420. mal_m[m] -> pyr[m]/ malic enzyme
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'mal-L[m]','pyr[m]'},[-1 -1; 0 1000]);
    model.c(ismember(model.rxns,'sink_pyr[m]'))=1;
    try
        eFBA=entropicFluxBalanceAnalysis(model,param);
        eFBA.f=sum(eFBA.v(find(model.c>0)));
    catch
        eFBA.f=NaN;
    end
    TestSolution(k,1) = eFBA.f;
    TestSolutionName{k,1} = 'mal_L[m] -> pyr[m]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 421. glu_L[c] -> gln_L[c]/ glutamine synthase
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'glu-L[c]','gln-L[c]'},[-1 -1; 0 1000]);
    model.c(ismember(model.rxns,'sink_gln-L[c]'))=1;
    try
        eFBA=entropicFluxBalanceAnalysis(model,param);
        eFBA.f=sum(eFBA.v(find(model.c>0)));
    catch
        eFBA.f=NaN;
    end
    TestSolution(k,1) = eFBA.f;
    TestSolutionName{k,1} = 'glu_L[c] -> gln_L[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 422. cys_L[c] -> coa[c]/ CoA synthesis from cysteine
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'cys-L[c]','coa[c]'},[-1 -1; 0 100]);
    model.lb(find(ismember(model.rxns,'sink_cys-L[c]')))=-1;
    model.ub(find(ismember(model.rxns,'sink_cys-L[c]')))=-1;
    model.lb(find(ismember(model.rxns,'sink_cys_L[c]')))=-1;
    model.ub(find(ismember(model.rxns,'sink_cys_L[c]')))=-1;
    model.c(ismember(model.rxns,'DPCOAK'))=1;
    try
        eFBA=entropicFluxBalanceAnalysis(model,param);
        eFBA.f=sum(eFBA.v(find(model.c>0)));
    catch
        eFBA.f=NaN;
    end
    if isempty(eFBA.f)
        eFBA.f=0;
    end
    TestSolution(k,1) = eFBA.f;
    TestSolutionName{k,1} = 'cys_L[c] -> coa[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 423. occoa[m] -> accoa[m]/ octanoate oxidation
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'occoa[m]','accoa[m]'},[-1 -1; 0 1000]);
    model.c(ismember(model.rxns,'sink_accoa[m]'))=1;
    try
        eFBA=entropicFluxBalanceAnalysis(model,param);
        eFBA.f=sum(eFBA.v(find(model.c>0)));
    catch
        eFBA.f=NaN;
    end
    TestSolution(k,1) = eFBA.f;
    TestSolutionName{k,1} = 'occoa[m] -> accoa[m]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 424. lnlncgcoa[c] -> dlnlcgcoa[c]/ fatty acid elongation
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'lnlncgcoa[c]','dlnlcgcoa[c]'},[-1 -1; 0 1000]);
    model.c(ismember(model.rxns,'sink_dlnlcgcoa[c]'))=1;
    try
        eFBA=entropicFluxBalanceAnalysis(model,param);
        eFBA.f=sum(eFBA.v(find(model.c>0)));
    catch
        eFBA.f=NaN;
    end
    TestSolution(k,1) = eFBA.f;
    TestSolutionName{k,1} = 'lnlncgcoa[c] -> dlnlcgcoa[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 425. chol[c] -> ach[c]/ acetyl-choline synthesis in brain
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'chol[c]','ach[c]'},[-1 -1; 0 1000]);
    model.c(ismember(model.rxns,'sink_ach[c]'))=1;
    try
        eFBA=entropicFluxBalanceAnalysis(model,param);
        eFBA.f=sum(eFBA.v(find(model.c>0)));
    catch
        eFBA.f=NaN;
    end
    TestSolution(k,1) = eFBA.f;
    TestSolutionName{k,1} = 'chol[c] -> ach[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 426. pyr[m] -> oaa[m]/ pyruvate carboxylase
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'pyr[m]','oaa[m]'},[-1 -1; 0 1000]);
    model.c(ismember(model.rxns,'sink_oaa[m]'))=1;
    try
        eFBA=entropicFluxBalanceAnalysis(model,param);
        eFBA.f=sum(eFBA.v(find(model.c>0)));
    catch
        eFBA.f=NaN;
    end
    TestSolution(k,1) = eFBA.f;
    TestSolutionName{k,1} = 'pyr[m] -> oaa[m]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 427. GABA aminotransferase
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_glu-L[e]'))=-1;model.ub(ismember(model.rxns,'EX_glu-L[e]'))=-1;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    model.c(ismember(model.rxns,'ABTArm'))=1;
    try
        eFBA=entropicFluxBalanceAnalysis(model,param);
        eFBA.f=sum(eFBA.v(find(model.c>0)));
    catch
        eFBA.f=NaN;
    end
    TestSolution(k,1) = eFBA.f;
    TestSolutionName{k,1} = 'GABA aminotransferase';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 428. methionine adenosyltransferase
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_met_L[e]'))=-1;model.ub(ismember(model.rxns,'EX_met_L[e]'))=-1;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    model.c(ismember(model.rxns,'METAT'))=1;
    try
        eFBA=entropicFluxBalanceAnalysis(model,param);
        eFBA.f=sum(eFBA.v(find(model.c>0)));
    catch
        eFBA.f=NaN;
    end
    TestSolution(k,1) = eFBA.f;
    TestSolutionName{k,1} = 'methionine adenosyltransferase';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 429. creatine synthesis
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_arg_L[e]'))=-1;model.ub(ismember(model.rxns,'EX_arg_L[e]'))=-1;
    model.lb(ismember(model.rxns,'EX_gly[e]'))=-1;model.ub(ismember(model.rxns,'EX_gly[e]'))=-1;
    [model] = addSinkReactions(model,{'crtn[c]'},[0 1000]);
    model.c(ismember(model.rxns,'sink_crtn[c]'))=1;
    try
        eFBA=entropicFluxBalanceAnalysis(model,param);
        eFBA.f=sum(eFBA.v(find(model.c>0)));
    catch
        eFBA.f=NaN;
    end
    TestSolution(k,1) = eFBA.f;
    TestSolutionName{k,1} = 'creatine synthesis';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 430. arachd[c] -> leuktrE4[c]/ leukotriene synthesis
    % requires multiple medium compounds --> RPMI
    
    model = modelOri;
    for i = 1 : length(RPMI_composition)
        model = changeRxnBounds(model,RPMI_composition{i},-1,'l');
    end
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'arachd[c]','leuktrE4[c]'},[-1 -1; 0 1000]);
    model.c(ismember(model.rxns,'sink_leuktrE4[c]'))=1;
    try
        eFBA=entropicFluxBalanceAnalysis(model,param);
        eFBA.f=sum(eFBA.v(find(model.c>0)));
    catch
        eFBA.f=NaN;
    end
    TestSolution(k,1) = eFBA.f;
    TestSolutionName{k,1} = 'arachd[c] -> leuktrE4[c] (with RPMI medium)';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 431. arachd[c] -> C06314[c]/ lipoxin synthesis
    model = modelOri;
    for i = 1 : length(RPMI_composition)
        model = changeRxnBounds(model,RPMI_composition{i},-1,'l');
    end
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'arachd[c]','C06314[c]'},[-1 -1; 0 1000]);
    model.c(ismember(model.rxns,'sink_C06314[c]'))=1;
    try
        eFBA=entropicFluxBalanceAnalysis(model,param);
        eFBA.f=sum(eFBA.v(find(model.c>0)));
    catch
        eFBA.f=NaN;
    end
    TestSolution(k,1) = eFBA.f;
    TestSolutionName{k,1} = 'arachd[c] -> C06314[c] (with RPMI medium)';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 432. nrpphr[c] -> 3mox4hoxm[c]/ degradation of norepinephrine
    model = modelOri;
    for i = 1 : length(RPMI_composition)
        model = changeRxnBounds(model,RPMI_composition{i},-1,'l');
    end
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'nrpphr[c]','3mox4hoxm[c]'},[-1 -1; 0 1000]);
    model.c(ismember(model.rxns,'sink_3mox4hoxm[c]'))=1;
    try
        eFBA=entropicFluxBalanceAnalysis(model,param);
        eFBA.f=sum(eFBA.v(find(model.c>0)));
    catch
        eFBA.f=NaN;
    end
    TestSolution(k,1) = eFBA.f;
    TestSolutionName{k,1} = 'nrpphr[c] -> 3mox4hoxm[c] (with RPMI medium)';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    %% 433. sbt_D[c] -> fru[c]/sorbitol pathway
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'sbt-D[c]','fru[c]'},[-1 -1; 0 1000]);
    model.c(ismember(model.rxns,'sink_fru[c]'))=1;
    try
        eFBA=entropicFluxBalanceAnalysis(model,param);
        eFBA.f=sum(eFBA.v(find(model.c>0)));
    catch
        eFBA.f=NaN;
    end
    TestSolution(k,1) = eFBA.f;
    TestSolutionName{k,1} = 'sbt_D[c] -> fru[c]/sorbitol pathway';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    
    %% 434. new addition 26.04.2017
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_glc[e]'))=-1;model.ub(ismember(model.rxns,'EX_glc[e]'))=-1;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    [model] = addSinkReactions(model,{'accoa[m]'},0,  1000);
    model.lb(find(ismember(model.rxns,'sink_coa[c]')))=-1;
    model.ub(find(ismember(model.rxns,'sink_coa[c]')))=1;
    model.c(ismember(model.rxns,'sink_accoa[m]'))=1;
    try
        eFBA=entropicFluxBalanceAnalysis(model,param);
        eFBA.f=sum(eFBA.v(find(model.c>0)));
    catch
        eFBA.f=NaN;
    end
    TestSolution(k,1) = eFBA.f;
    TestSolutionName{k,1} = 'Mitochondrial accoa de novo synthesis from glc';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
    
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_glc[e]'))=-1;model.ub(ismember(model.rxns,'EX_glc[e]'))=-1;
    model.lb(ismember(model.rxns,'EX_o2[e]'))=-40;model.ub(ismember(model.rxns,'EX_o2[e]'))=-1;
    [model] = addSinkReactions(model,{'succoa[m]'},0,1000);
    model.c(ismember(model.rxns,'sink_succoa[m]'))=1;
    model.lb(find(ismember(model.rxns,'sink_coa[c]')))=-1;
    model.ub(find(ismember(model.rxns,'sink_coa[c]')))=1;
    try
        eFBA=entropicFluxBalanceAnalysis(model,param);
        eFBA.f=sum(eFBA.v(find(model.c>0)));
    catch
        eFBA.f=NaN;
    end
    TestSolution(k,1) = eFBA.f;
    TestSolutionName{k,1} = 'Mitochondrial succoa de novo synthesis from glc';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(eFBA.v)>tol))]; end ;k = k +1;clear eFBA
end
TestSolution(find(abs(TestSolution)<tol))=0;
TestSolutionName(:,2) = num2cell(TestSolution);
TestedRxns = unique(TestedRxns);
TestedRxns = intersect(modelOri.rxns,TestedRxns); % only those reactions that are also in modelOri not those that have been added to the network
PercTestedRxns = length(TestedRxns)*100/length(modelOri.rxns);
if saveDiary
    diary off
end
