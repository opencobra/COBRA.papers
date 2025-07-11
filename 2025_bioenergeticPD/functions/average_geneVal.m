C1 = readtable('expressionFile_candidateSamples_Control_1');
C2 = readtable('expressionFile_candidateSamples_Control_2');
P1 = readtable('expressionFile_candidateSamples_PD_1');
P2 = readtable('expressionFile_candidateSamples_PD_2');

% Import the gene libraries and adjust the length for NaN
addpath('~/work/sbgCloud/programReconstruction/projects/bioenergetics/results/modelGeneration/inputs/omics/transcriptomics/datasets_GEO/processedOmics/xlsx/');
hgu133a = xlsread('library_Probe_Symbol_Entrez_IDs.xlsx', 'HG_U133A');
hgu133a = [NaN; hgu133a];
if length(hgu133a) < length(C1.Var1);
    hgu133a(end:length(C1.Var1)) = [NaN];
end
hgu133plus2 = xlsread('library_Probe_Symbol_Entrez_IDs.xlsx', 'HG_U133Plus2');
hgu133plus2 = [NaN; hgu133plus2];
if length(hgu133plus2) < length(C2.Var1);
    hgu133plus2(end:length(C2.Var1)) = [NaN];
end

% Add probe IDs
[num, probeshgu133a] = xlsread('library_Probe_Symbol_Entrez_IDs.xlsx',...
    'HG_U133A', 'A2:A22284');
[num, probeshgu133plus2] = xlsread('library_Probe_Symbol_Entrez_IDs.xlsx',...
    'HG_U133Plus2', 'A2:A54676');
C1.Var1 = probeshgu133a;
C2.Var1 = probeshgu133plus2;
P1.Var1 = probeshgu133a;
P2.Var1 = probeshgu133plus2;

% Replace Entrez IDs for quality control
C1.Var2 = hgu133a;
C2.Var2 = hgu133plus2;
P1.Var2 = hgu133a;
P2.Var2 = hgu133plus2;

% Remove probes with undesirable probe set extensions
check = {'AFFX-','_x_at', '_s_at', '_a_at', '_r_at', '_st', ...
    '_i', '_g_at', '_f_at', '_i_at', '_b_at', '_l_at', '_r_at'};

for i = 1:length(C1.Var1);
    probeFound = ~cellfun(@isempty,regexp(C1.Var1{i},check));
    if find(probeFound)>0;
        C1.Var3(i) = 0; %remove
    else C1.Var3(i) = 1; %keep
    end
end

for i = 1:length(C2.Var1);
    probeFound = ~cellfun(@isempty,regexp(C2.Var1{i},check));
    if find(probeFound)>0;
        C2.Var3(i) = 0; %remove
    else C2.Var3(i) = 1; %keep
    end
end

for i = 1:length(P1.Var1);
    probeFound = ~cellfun(@isempty,regexp(P1.Var1{i},check));
    if find(probeFound)>0;
        P1.Var3(i) = 0; %remove
    else P1.Var3(i) = 1; %keep
    end
end

for i = 1:length(P2.Var1);
    probeFound = ~cellfun(@isempty,regexp(P2.Var1{i},check));
    if find(probeFound)>0;
        P2.Var3(i) = 0; %remove
    else P2.Var3(i) = 1; %keep
    end
end


% Remove unwanted probe IDs:
C1 = sortrows(C1,'Var3','ascend');
C2 = sortrows(C2,'Var3','ascend');
P1 = sortrows(P1,'Var3','ascend');
P2 = sortrows(P2,'Var3','ascend');

C1(1:10509,:) = [];
C2(1:17320,:) = [];
P1(1:10509,:) = [];
P2(1:17320,:) = [];

% Move Entrez IDs next to probe IDs
C1 = C1(:,[1 end-1 2:end]);
C2 = C2(:,[1 end-1 2:end]);
P1 = P1(:,[1 end-1 2:end]);
P2 = P2(:,[1 end-1 2:end]);

% Delete variables that are no longer required
C1(:,{'Var1','Var3','Var2_1'}) = [];
C2(:,{'Var1','Var3','Var2_1'}) = [];
P1(:,{'Var1','Var3','Var2_1'}) = [];
P2(:,{'Var1','Var3','Var2_1'}) = [];

% Remove the rows that were not mapped to an Entrez ID
C1 = sortrows(C1,'Var2','ascend');
C2 = sortrows(C2,'Var2','ascend');
P1 = sortrows(P1,'Var2','ascend');
P2 = sortrows(P2,'Var2','ascend');

C1(10429:end,:) = [];
C2(27414:end,:) = [];
P1(10429:end,:) = [];
P2(27414:end,:) = [];

% % remove NaN rows
% C1(19855:end,:) = [];
% C2(42323:end,:) = [];
% P1(19855:end,:) = [];
% P2(42323:end,:) = [];

% Remove samples with a GNUSE > 1.25 (for HGU133a & HGU133plus2) as per
% McCall et al., Assessing affymetrix GeneChip microarray quality., BMC
% Bioinformatics. 2011
C1(:,{'GSM508692_CEL_gz','GSM506020_CEL_gz','GSM508703_CEL_gz',...
    'GSM508723_CEL_gz','GSM508734_CEL_gz'}) = [];
C2(:,{'GSM503950_C_1074p_SNc_CEL_gz','GSM503951_C_1271p_SNc_CEL_gz',...
    'GSM503952_C_2829_SNc_CEL_gz','GSM503953_C_3132p_SNc_CEL_gz',...
    'GSM503954_C_3397_SNc_CEL_gz','GSM503955_C_3543_SNc_CEL_gz',...
    'GSM503956_C_3603_SNc_CEL_gz','GSM503957_C_5220_SNc_CEL_gz'}) = [];
P1(:,{'GSM508628_CEL_gz','GSM508732_CEL_gz'}) = [];
P2(:,{'GSM503958_PD_1364_SNc_CEL_gz','GSM503959_PD_1401_SNc_CEL_gz',...
    'GSM503962_PD_2525_SNc_CEL_gz','GSM503963_PD_3769p_SNc_CEL_gz',...
    'GSM503964_PD_3790p_SNc_CEL_gz','GSM503960_PD_1647_SNc_CEL_gz',...
    'GSM503961_PD_2515p_SNc_CEL_gz','GSM503966_PD_5138_SNc_CEL_gz',...
    'GSM503967_PD_5476_SNc_CEL_gz','GSM503965_PD_3803p_SNc_CEL_gz'}) = [];

% Obtain unique Entrez IDs and their corresponding gene ExVal. Take the mean of the gene Exval

% Extract the unique genes and symbols for platforms hgu133a and
% hgu133plus2
genes1 = unique(C1{:,1});
genes2 = unique(C2{:,1});

% Merge all the probes by taking the mean of the corresponding gene ExVal.
C1tmp = genes1;
samplesC1 = {};
for j = 3:size(C1,2);
    for i = 1:length(C1tmp);
        tmp = find(ismember(C1{:,1},C1tmp(i)));
        if size(tmp,1) == 1;
            C1tmp(i,j-1) = C1{tmp,j};
        elseif size(tmp,1) >1;
            C1tmp(i,j-1) = mean(C1{tmp,j});
        end
    end
    samplesC1{j-2} = C1.Properties.VariableNames{j};
end

C2tmp = genes2;
samplesC2 = {};
for j = 3:size(C2,2);
    for i = 1:length(C2tmp);
        tmp = find(ismember(C2{:,1},C2tmp(i)));
        if size(tmp,1) == 1;
            C2tmp(i,j-1) = C2{tmp,j};
        elseif size(tmp,1) >1;
            C2tmp(i,j-1) = mean(C2{tmp,j});
        end
    end
    samplesC2{j-2} = C2.Properties.VariableNames{j};
end

P1tmp = genes1;
samplesP1 = {};
for j = 3:size(P1,2);
    for i = 1:length(P1tmp);
        tmp = find(ismember(P1{:,1},P1tmp(i)));
        if size(tmp,1) == 1;
            P1tmp(i,j-1) = P1{tmp,j};
        elseif size(tmp,1) >1;
            P1tmp(i,j-1) = mean(P1{tmp,j});
        end
    end
    samplesP1{j-2} = P1.Properties.VariableNames{j};
end

P2tmp = genes2;
samplesP2 = {};
for j = 3:size(P2,2);
    for i = 1:length(P2tmp);
        tmp = find(ismember(P2{:,1},P2tmp(i)));
        if size(tmp,1) == 1;
            P2tmp(i,j-1) = P2{tmp,j};
        elseif size(tmp,1) >1;
            P2tmp(i,j-1) = mean(P2{tmp,j});
        end
    end
    samplesP2{j-2} = P2.Properties.VariableNames{j};
end

% Obtain the average gene ExVal of each row (i.e. of each Entrez Gene ID)
averageC1 = C1tmp(:,1);
j = 0;
for i = 1:length(C1tmp);
    j= j +1;
    averageC1(j,2) = mean(C1tmp(i,2:end));
end

averageC2 = C2tmp(:,1);
j = 0;
for i = 1:length(C2tmp);
    j= j +1;
    averageC2(j,2) = mean(C2tmp(i,2:end));
end

averageP1 = P1tmp(:,1);
j = 0;
for i = 1:length(P1tmp);
    j= j +1;
    averageP1(j,2) = mean(P1tmp(i,2:end));
end

averageP2 = P2tmp(:,1);
j = 0;
for i = 1:length(P2tmp);
    j= j +1;
    averageP2(j,2) = mean(P2tmp(i,2:end));
end

% Now combine the two control and the two PD sets by taking once again the
% average (i.e. average of C1 & C2 and the average of P1 & P2)
allGenesC = unique([averageC1(:,1); averageC2(:,1)]);
allGenesP = unique([averageP1(:,1); averageP2(:,1)]);

summaryC = allGenesC;
for i = 1:length(summaryC);
    k = summaryC(i,1);
    tmp1 = find(ismember(averageC1(:,1), k));
    tmp2 = find(ismember(averageC2(:,1), k));
    if size(tmp1,2) == 1 && size(tmp2,2) == 1;
        summaryC(i,2) = mean([averageC1(tmp1,2), averageC2(tmp2,2)]);
    elseif size(tmp1,2) == 1;
        summaryC(i,2) = averageC1(tmp1,2);
    elseif size(tmp2,2) == 1;
        summaryC(i,2) = averageC2(tmp2,2);
    end
end
summaryC = table(summaryC(:,1), summaryC(:,2), 'VariableNames', {'ID','values'});
summaryC = sortrows(summaryC,'values','descend');

summaryP = allGenesP;
for i = 1:length(summaryP);
    k = summaryP(i,1);
    tmp1 = find(ismember(averageP1(:,1), k));
    tmp2 = find(ismember(averageP2(:,1), k));
    if size(tmp1,2) == 1 && size(tmp2,2) == 1;
        summaryP(i,2) = mean([averageP1(tmp1,2), averageP2(tmp2,2)]);
    elseif size(tmp1,2) == 1;
        summaryP(i,2) = averageP1(tmp1,2);
    elseif size(tmp2,2) == 1;
        summaryP(i,2) = averageP2(tmp2,2);
    end
end
summaryP = table(summaryP(:,1), summaryP(:,2), 'VariableNames', {'ID','values'});
summaryP = sortrows(summaryP,'values','descend');