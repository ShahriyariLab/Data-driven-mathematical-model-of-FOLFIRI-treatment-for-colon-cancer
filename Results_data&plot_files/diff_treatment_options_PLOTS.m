%Different treatment dynamics 
% diff_treatment_options_PLOTS: script for plotting the results in Figure 4 obtained by running 
%   diff_treatment_options.ipynb

% Author: Aparajita Budithi, Arkadz Kirshtein, https://sites.google.com/site/akirshtein/
% (c) Shahriyari Lab https://sites.google.com/site/leilishahriyari/
%
% If using this or related code please cite 
% Budithi, A.; Su, S.; Kirshtein, A.; Shahriyari L. 
%   Data driven mathematical model of FOLFIRI treatment for colon cancer. Cancers2021. 
%   (Manuscript submitted for publishing)


clear;
% total time = start-treatment + treatment + observation time 
%T= 0 + 84 + 365 ;
vars={'Naive T-cells', 'helper T-cells', 'cytotoxic cells', 'Treg-cells', 'Naive Dendritic cells', 'Dendritic cells', 'Macrophages', 'Cancer', 'Necrotic cells', '\mu_1', '\mu_2', 'HMGB1', 'IFN-\gamma', 'TGF-\beta', '5FU', 'IR', 'LV', 'Total cells'};
tab=readtable('input/Large_Tumor_variable_scaling.csv');
cells=[tab{:,:}];
lmod=1:9; %all cells
T = 1000;
data=cell(1,5);
drug_comb = {'FU','FULV','FUIR','FUIRLV'};
for k =1:4
    for cluster=1:5
        tab=readtable(['Results/Diff_treatment/' drug_comb{k} '/V63-cluster-' num2str(cluster) '-of-5-results-dat.csv']);
        data{cluster}=tab{:,:};
    end
    clear cluster t tab
    save(['Results\Diff_treatment\' drug_comb{k} '\V63-results-all-data.mat'])
% reading parsed data
    load(['Results\Diff_treatment\' drug_comb{k} '\V63-results-all-data.mat'])
% total time = start-treatment + treatment + observation time 
    figure;
    
    cmap=lines(5);
    for i=1:9
        subplot(2,5,i);
        hold on
        for cluster=1:5
            plot(data{cluster}(:,1),data{cluster}(:,i+1)*cells(cluster,i),'color',cmap(cluster,:));
        end
        xlabel('time (days)')
        ylabel(vars{i})
        xlim([0 T])
        hold off
    end
    %total cells
    subplot(2,5,10);
    hold on
    for cluster=1:5
        plot(data{cluster}(:,1),data{cluster}(:,2:10)*cells(cluster,lmod)','color',cmap(cluster,:));
    end
    xlabel('time (days)')
    ylabel(vars{end})
    xlim([0 T])
    hold off
    legend({'cluster 1','cluster 2','cluster 3','cluster 4','cluster 5'})
    sgtitle(drug_comb{k})
end
