% Varying Start-treatment dynamics - Figure 5
% Data file link:
https://umass-my.sharepoint.com/:f:/g/personal/abudithi_umass_edu/EiCO4DSSWG5Cj_5plwdvs3QBml6hpPKF85AIBc447AdRKQ?e=a4J7Mo


clear;
% total time = start-treatment + treatment + observation time 
%T= 0 + 84 + 365 ;
vars={'Naive T-cells', 'helper T-cells', 'cytotoxic cells', 'Treg-cells', 'Naive Dendritic cells', 'Dendritic cells', 'Macrophages', 'Cancer', 'Necrotic cells', '\mu_1', '\mu_2', 'HMGB1', 'IFN-\gamma', 'TGF-\beta', '5FU', 'IR', 'LV', 'Total cells'};
tab=readtable('input/Large_Tumor_variable_scaling.csv');
cells=[tab{:,:}];
lmod=1:9; %all cells
% total time  
T=5000;
data=cell(1,5);
start_treatment_t = [365, 365*3, 365*5, 365*7] ;
for k =1:length(start_treatment_t)
    time = start_treatment_t(k);
    for cluster=1:5
        tab=readtable(['Results/Var_start_treatment/' num2str(time) '-days-cluster-' num2str(cluster) '-of-5-results-dat.csv']);
        data{cluster}=tab{:,:};
    end
    save(['Results\Var_start_treatment\' num2str(time) '-days-results-all-data.mat'])
    
    clear cluster t tab

    % reading parsed data
    load(['Results\Var_start_treatment\' num2str(time) '-days-results-all-data.mat'])
    
    figure;
    cmap=lines(5);
    for i=1:9
        subplot(2,5,i);
        hold on
        for cluster=1:5
            plot(data{cluster}(:,1),data{cluster}(:,i+1)*cells(cluster,i),'color',cmap(cluster,:),'LineWidth',1);
        end
        xlabel('time (days)','FontSize',12)
        ylabel(vars{i},'FontSize',12)
        xlim([0 T])
        hold off
    end
    %total cells
    subplot(2,5,10);
    hold on
    for cluster=1:5
        plot(data{cluster}(:,1),data{cluster}(:,2:10)*cells(cluster,lmod)','color',cmap(cluster,:),'LineWidth',1);
    end
    xlabel('time (days)','FontSize',12)
    ylabel(vars{end},'FontSize',12)
    xlim([0 T])
    hold off
    legend({'cluster 1','cluster 2','cluster 3','cluster 4','cluster 5'})
    sgtitle(['Start treatment- ' num2str(time/365) ' years'])
    
end
