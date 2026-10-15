function renamedVmhList = renameVmhToMetName_2026Gout(vmhList, reverse)

if nargin <2
    reverse = false;
end

tPose = false;
if size(vmhList,1)==1 && numel(vmhList)~=1
    vmhList = vmhList';
    tPose = true;
end

% Process VMH metabolites
% vmhList = string(vmhList);
vmhMets = erase(vmhList,["DM_","[bc]","_bc_"]);

% Load VMH database
metabolites = loadVMHDatabase().metabolites;

% Find metabolite names
[~,ia,ib] = intersect(vmhMets,metabolites(:,1),'stable');
metaboliteList = table(vmhList(ia),metabolites(ib,2),'VariableNames',{'VMH','Metabolite'});

% 
renamedVmhList = vmhList;
for i=1:numel(metaboliteList.VMH)
    renamedVmhList = replace(renamedVmhList, metaboliteList.VMH(i), metaboliteList.Metabolite(i));
end

renamedVmhList = replace(renamedVmhList, 'Deoxycholic acid-24glucuronide, CDA-24G', 'Deoxycholic acid-24G');
renamedVmhList = replace(renamedVmhList, '3,4-Dihydroxyphenylethyleneglycol', '3,4-Dihydroxyphenylglycol');
renamedVmhList = replace(renamedVmhList, 'creatine', 'Creatine');
renamedVmhList = replace(renamedVmhList, 'gDeoxycholic acid 3-sulfate', 'Glycodeoxycholic acid 3-sulfate');
renamedVmhList = replace(renamedVmhList, 'deoxycholic acid', 'Deoxycholic acid');
renamedVmhList = replace(renamedVmhList, 'TauroDeoxycholic acid 3-sulfate', 'Taurodeoxycholic acid 3-sulfate');
renamedVmhList = replace(renamedVmhList, 'taurolithocholate', 'Taurolithocholate');
renamedVmhList = replace(renamedVmhList, 'GlycoDeoxycholic acid 3-sulfate', 'Glycodeoxycholic acid 3-sulfate');
renamedVmhList = replace(renamedVmhList, 'formaldehyde', 'Formaldehyde');





if reverse == true
    % Rename metabolites to VMH reactions
    vmhList = string(vmhList);
    vmhList(matches(vmhList,"S-Adenosyl-L-methionine")) = "DM_amet[bc]";
    vmhList(matches(vmhList,"L-arginine")) = "DM_arg_L[bc]";
    vmhList(matches(vmhList,"Creatine")) = "DM_creat[bc]";
    vmhList(matches(vmhList,"Taurine")) = "DM_taur[bc]";
    vmhList(matches(vmhList,"L-lysine")) = "DM_lys_L[bc]";
    vmhList(matches(vmhList,"Deoxycholic acid-24G")) = "DM_dca24g[bc]";
    vmhList(matches(vmhList,"Taurodeoxycholic acid 3-sulfate")) = "DM_tdca3s[bc]";
    vmhList(matches(vmhList,"formaldehyde")) = "DM_fald[bc]";
    vmhList(matches(vmhList,"Formate")) = "DM_for[bc]";
    vmhList(matches(vmhList,"deoxycholic acid")) = "DM_dchac[bc]";
    vmhList(matches(vmhList,"3-dehydro-Deoxycholate")) = "DM_3dhchol[bc]";
    vmhList(matches(vmhList,"3-Dehydrocholate")) = "DM_3dhchol[bc]";
    vmhList(matches(vmhList,"ChenoDeoxycholic acid-3glucuronide, CDCA-3G")) = "DM_cdca3g[bc]";
    vmhList(matches(vmhList,"Deoxycholic acid 3-sulfate")) = "DM_dca3s[bc]";
    vmhList(matches(vmhList,"Glycodeoxycholic acid 3-sulfate")) = "DM_gdca3s[bc]";
    vmhList(matches(vmhList,"D-glucose")) = "DM_glc_D[bc]";
    vmhList(matches(vmhList,"Isocholate")) = "DM_isochol[bc]";
    vmhList(matches(vmhList,"Propionate")) = "DM_ppa[bc]";
    vmhList(matches(vmhList,"Formaldehyde")) = "DM_fald[bc]";
    
    
    renamedVmhList = cellstr(vmhList);
end


if tPose == true
    renamedVmhList = renamedVmhList';
end

% % Rename VMH reactions to metabolites
% vmhList = string(vmhList);
% vmhList(matches(vmhList,"DM_leuleu[bc]")) = "Leucylleucine";
% vmhList(matches(vmhList,"DM_leu_L[bc]")) = "L-leucine";
% vmhList(matches(vmhList,"DM_but[bc]")) = "Butyrate";
% vmhList(matches(vmhList,"DM_pnto_R[bc]")) = "Pantothenate";
% vmhList(matches(vmhList,"DM_nac[bc]")) = "Nicotinic acid";
% vmhList(matches(vmhList,"DM_ttdca[bc]")) = "Myristic acid";
% metaboliteList = cellstr(vmhList);
% 
% if reverse == true
%     % Rename metabolites to VMH reactions
%     vmhList = string(vmhList);
%     vmhList(matches(vmhList,"Leucylleucine")) = "DM_leuleu[bc]";
%     vmhList(matches(vmhList,"L-leucine")) = "DM_leu_L[bc]";
%     vmhList(matches(vmhList,"Butyrate")) = "DM_but[bc]";
%     vmhList(matches(vmhList,"Pantothenate")) = "DM_pnto_R[bc]";
%     vmhList(matches(vmhList,"Nicotinic acid")) = "DM_nac[bc]";
%     vmhList(matches(vmhList,"Myristic acid")) = "DM_ttdca[bc]";
%     metaboliteList = cellstr(vmhList);
% 
% end