%choose the input data
cd(modelCollectionDirectory)

load('~/work/sbgCloud/data/models/published/Recon1_reconstruction/09.17.09.Recon1.mat');
Recon1 = modelSBML;clear modelSBML;
Recon1.modelID='Recon1.0';
save('Recon1.0.mat','Recon1');

if 0
load('~/work/sbgCloud/data/models/published/Recon2_reconstruction/121114_Recon2betaRecon.mat')
Recon2 = modelRecon2beta121114Recon; clear modelRecon2beta121114Recon;
Recon2.modelID='Recon2.0';
save('Recon2.0.mat','Recon2');

load('~/work/sbgCloud/data/models/published/Recon2_model/121114_Recon2betaModel.mat')
Recon2model = modelRecon2beta121114; clear modelRecon2beta121114;
Recon2model.modelID='Recon2.0model';
save('Recon2.0model.mat','Recon2model');


load('~/work/sbgCloud/data/models/published/Recon2.02/recon2.v02.mat')
Recon202 = model; clear model;
Recon202.modelID='Recon2.02';
save('Recon2.02.mat','Recon202');

load('~/work/sbgCloud/data/models/published/Recon2.02/recon2model.v02.mat')
Recon202model = model; clear model;
Recon202model.modelID='Recon2.02model';
save('Recon2.02model.mat','Recon202model');
end

load('~/work/sbgCloud/data/models/published/Recon2.04_model/Recon2.v04.mat');
Recon204model=model;clear model;
Recon204model.modelID='Recon2.04model';
save('Recon2.04model.mat','Recon204model')

load('~/work/sbgCloud/data/models/published/Recon2.2/Recon2_2_biomodel.mat');
Recon22 = model; clear model;
Recon22.modelID='Recon2.2model';
save('Recon2.2model.mat','Recon22')

load('~/work/sbgCloud/data/models/published/HMR2/HMRdatabase2_00.mat');
HMR2 = ihuman; clear ihuman;
HMR2.modelID='HMR2.0';
HMR2 = findSExRxnInd(HMR2);
bool=endsWith(HMR2.mets,'x');
metaboliteList=HMR2.mets(bool);
removeRxnFlag=0;
HMR2 = removeMetabolites(HMR2,metaboliteList,removeRxnFlag);
save('HMR2.0.mat','HMR2')

load('~/work/sbgCloud/data/models/published/Recon3D_301_Thiele_2018/Recon3D_301.mat');
Recon301 = Recon3D;
Recon3.modelID='Recon3.01';
save('Recon3.01.mat','Recon301')

load('~/work/sbgCloud/data/models/published/Recon3D_301_Thiele_2018/Recon3DModel_301.mat');
Recon301model = Recon3DModel;
Recon301model.modelID='Recon3.01model';
save('Recon3.01model.mat','Recon301model')

%humanGEM1 = readCbModel('~/work/sbgCloud/code/Human-GEM-1.0.0/ModelFiles/xml/humanGEM.xml');
%The model contains 68 errors and 1 warnings.
%Error encountered during read.
% https://github.com/SysBioChalmers/RAVEN/blob/main/core/getExchangeRxns.m
load('~/work/sbgCloud/data/models/published/Human-GEM-1.0.0/ModelFiles/mat/humanGEM1.mat');
humanGEM1.modelID='humanGEM1';
humanGEM1.S = sparse(humanGEM1.S);
save('humanGEM1.mat','humanGEM1')

if 0
%sha = 26005b7de4abcce9edd7dcf3330c2075b36245fb
%humanGEM1p10 = readCbModel('/home/rfleming/work/sbgCloud/code/Human-GEM/model/Human-GEM.xml');
%The model contains 0 errors and 1 warnings.
load('~/work/sbgCloud/data/models/published/Human-GEM/model/humanGEM1p10.mat');
humanGEM1p10.modelID='humanGEM1p10';
humanGEM1p10.S = sparse(humanGEM1p10.S);
save('humanGEM1p10.mat','humanGEM1p10')
end

if 0
load('~/work/sbgCloud/data/models/published/Harvetta/Harvetta_1_03c.mat')
Harvetta = female;
save('Harvetta.mat','Harvetta');

load('~/work/sbgCloud/data/models/published/Harvey/Harvey_1_03c.mat')
Harvey = male;
save('Harvey.mat','Harvey');
end
