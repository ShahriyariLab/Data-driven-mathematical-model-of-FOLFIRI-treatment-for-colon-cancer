
% var_sensitive_par_PLOTS: script for plotting the results in Figure 9 & 10, obtained by running 
%   folfiri_dynamics.ipynb and ParseQSPData.m
% Author: Aparajita Budithi
% (c) Shahriyari Lab https://sites.google.com/site/leilishahriyari/
%
% If using this or related code please cite 
% Budithi, A.; Su, S.; Kirshtein, A.; Shahriyari L. 
%   Data driven mathematical model of FOLFIRI treatment for colon cancer. Cancers2021. 
%   (Manuscript submitted for publishing)

clear;
% Scaling
vars={'Naive T-cells', 'helper T-cells', 'cytotoxic cells', 'Treg-cells', 'Naive Dendritic cells', 'Dendritic cells', 'Macrophages', 'Cancer', 'Necrotic cells', '\mu_1', '\mu_2', 'HMGB1', 'IFN-\gamma', 'TGF-\beta', '5FU', 'IR', 'LV', 'Total cells'};
tab=readtable('input/Large_Tumor_variable_scaling.csv');
cells=[tab{:,:}];
lmod=1:9; %all cells

params = {'delta_{CFULV}','delta_{5fuM}'};
param_values = [0.0005362,3.6166e-06];
% delta_5fuM = 3.6166e-06, delta_C5fuLV = 0.0005362

cons = [0.1,0.5,1,2,3,4,5];
t_points = [100, 169, 365, 365*2, 365*3];

%color variation for the plots
col = linspace(1,0,length(t_points));
col_2 = linspace(0,1,length(t_points)); 


for par = 1:2
    param = params{par};
    original_par = param_values(par);
    
    %initiate arrays for plot
    x_par = cons*original_par;
    Can_values = zeros(1,length(cons));
    Tot_cells = zeros(1,length(cons));
    Cancer_by_t_point = zeros(length(t_points),length(cons));
    Tot_cells_by_t_point = zeros(length(t_points),length(cons));
    
    for cluster = 1:5
        for k = 1:length(t_points)
            t = t_points(k);
            for i=1:length(cons)
                data = readtable(['Results/Var_sens_par/' param '/cluster-' num2str(cluster) '-results-' num2str(i-1) '-dat.csv']);
                pos = find(data.time == t);
                %pos = pos(1);
                Can_values(i) = data.Cancer(pos)*cells(cluster,8);
                Tot_cells(i) = table2array(data(pos,2:10))*transpose(cells(cluster,lmod));
            end
            Cancer_by_t_point(k,:) = Can_values;
            Tot_cells_by_t_point(k,:) = Tot_cells;
        end
        %PLOTS
        %index for original parameter value
        index = find(cons == 1);
        %Cancer
        figure;
        cmap=lines(length(t_points));
        for m = 1:length(t_points)
            hold on
            p(m) = plot(x_par,Cancer_by_t_point(m,:),'-o','color',[0,col(m),col_2(m)],'LineWidth',1.5,'MarkerSize',12);
            q = plot(x_par(index),Cancer_by_t_point(m,index),'k*','MarkerSize',12);
            hold off
            xlabel(param,'FontSize',18)
            ylabel(vars{8},'FontSize',18)
        end
        title(['Cluster ' num2str(cluster)],'FontSize',18)
        leg = legend([p,q],'100 days','169 days (end treatment)','1 year',...
            '2 years','3 years','original parameter value');
        title(leg,'Time point');
        %set(gca,'FontSize',18)

        %Total cells
        figure;
        cmap=lines(length(t_points));
        for m = 1:length(t_points)
            hold on
            p(m) = plot(x_par,Tot_cells_by_t_point(m,:),'-o','color',[0,col(m),col_2(m)],'LineWidth',1.5,'MarkerSize',12);
            q = plot(x_par(index),Tot_cells_by_t_point(m,index),'k*','MarkerSize',12);
            hold off
            xlabel(param,'FontSize',18)
            ylabel(vars{end},'FontSize',18)
        end
        title(['Cluster ' num2str(cluster)],'FontSize',18)
        leg = legend([p,q],'100 days','169 days (end treatment)','1 year',...
            '2 years','3 years','original parameter value');
        title(leg,'Time point');
        %set(gca,'FontSize',18)

    end
end
