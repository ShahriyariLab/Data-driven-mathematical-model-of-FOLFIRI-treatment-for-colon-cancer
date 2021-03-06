% ParseQSPData: script to parse the results obtained by running 
%   folfiri_dynamics.ipybn for plotting by PlotQSPData.m
% Author: Arkadz Kirshtein, https://sites.google.com/site/akirshtein/
% (c) Shahriyari Lab https://sites.google.com/site/leilishahriyari/
%
% part of https://github.com/ShahriyariLab/Data-driven-mathematical-model-of-FOLFIRI-treatment-for-colon-cancer
% If using this or related code please cite 
% Budithi, A.; Su, S.; Kirshtein, A.; Shahriyari L. 
%   Data driven mathematical model of FOLFIRI treatment for colon cancer. Cancers2021. 
%(Manuscript submitted for publishing)


%Parsing dynamics data

clear;
% total time = start-treatment + treatment + observation time 
%T= 0 + 84 + 365 ;
vars={'Naive T-cells', 'helper T-cells', 'cytotoxic cells', 'Treg-cells', 'Naive Dendritic cells', 'Dendritic cells', 'Macrophages', 'Cancer', 'Necrotic cells', '\mu_1', '\mu_2', 'HMGB1', 'IFN-\gamma', 'TGF-\beta', '5FU', 'IR', 'LV', 'Total cells'};
tab=readtable('input/Large_Tumor_variable_scaling.csv');
cells=[tab{:,:}];
lmod=1:9; %all cells

data=cell(1,5);
% maxdata=cell(1,5);
% mindata=cell(1,5);
for cluster=1:5
    tab=readtable(['Results/All_Dynamic/V63-cluster-' num2str(cluster) '-of-5-results-dat.csv']);
    data{cluster}=tab{:,:};
%     tab=readtable(['Data/Dynamic/V63-cluster-' num2str(cluster) '-of-5-results-varmindat.csv']);
%     mindata{cluster}=tab{:,:};
%     tab=readtable(['Data/Dynamic/V63-cluster-' num2str(cluster) '-of-5-results-varmaxdat.csv']);
%     maxdata{cluster}=tab{:,:};
end

clear cluster t tab

save('Results\All_Dynamic\V63-results-all-data.mat')
