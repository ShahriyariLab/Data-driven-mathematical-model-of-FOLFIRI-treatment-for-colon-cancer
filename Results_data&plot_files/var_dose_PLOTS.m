% var_dose_PLOTS: script for plotting the results in Figure 7 obtained by running 
%   var_dose.ipynb
% Author: Aparajita Budithi
% (c) Shahriyari Lab https://sites.google.com/site/leilishahriyari/
%
% If using this or related code please cite 
% Budithi, A.; Su, S.; Kirshtein, A.; Shahriyari L. 
%   Data driven mathematical model of FOLFIRI treatment for colon cancer. Cancers2021. 
%   (Manuscript submitted for publishing)
%Cancer & Total cells vs varying dosages

clear;

% Scaling
vars={'Naive T-cells', 'helper T-cells', 'cytotoxic cells', 'Treg-cells', 'Naive Dendritic cells', 'Dendritic cells', 'Macrophages', 'Cancer', 'Necrotic cells', '\mu_1', '\mu_2', 'HMGB1', 'IFN-\gamma', 'TGF-\beta', '5FU', 'IR', 'LV', 'Total cells'};
tab=readtable('input/Large_Tumor_variable_scaling.csv');
cells=[tab{:,:}];
lmod=1:9; %all cells

id = {'TCGA-CM-6172','TCGA-A6-6141','TCGA-A6-6142','TCGA-G4-6303','TCGA-A6-6781'};

%figure;
for j = 1:5
    %read patient data 
    pat_data = readtable('input/data_with_FIRST_end_day.csv');
    last_data = readtable('input/data_with_end_day.csv');
    num = find(string(pat_data.Patient_ID) == id{j});
    tumor_status = pat_data.tumor_status(num);
    first_follow_up = pat_data.end_day(num);
    last_follow_up = last_data.end_day(num);
    cluster = pat_data.Cluster(num);    
    drugs = {'LV','FU'};

    for k = 1:2
        
        if k == 1
            treatment_time = pat_data.LV_days_to_drug_therapy_end(num) - pat_data.LV_days_to_drug_therapy_start(num);
            cycle_time = treatment_time/pat_data.LV_number_cycles(num);
            prescribed_dose = pat_data.LV_initial(num);
            drug_name = 'Leucovorin';

        else 
            treatment_time = pat_data.FU_days_to_drug_therapy_end(num) - pat_data.FU_days_to_drug_therapy_start(num);
            cycle_time = treatment_time/pat_data.FU_number_cycles(num);
            prescribed_dose = pat_data.FU_initial(num);
            drug_name = 'Fluorouracil';

        end
        drug = drugs{k};
        %%PLOTS%%%

        cons = [0.1,0.5,1,2,5,7,10];
        drug_values = cons(:).*prescribed_dose;
        drug_values = round(drug_values,2);
        lmod = 1:9;
        col = linspace(0.7,0,length(cons));
        if id{j} == 'TCGA-G4-6303'
            T = 6500;
        else
            T = 3500;
        end
        
        %Cancer
        figure;
        for i = 1:length(cons)
            data=readtable(['Results/Var_doses/Patient_' id{j}(9:12) '_' tumor_status{1} '/' drug '-' num2str(i-1) '-dat.csv']);
            pos = find(data.time == last_follow_up);
            pos1 = find(data.time == first_follow_up);
            hold on
            if i == 3
                p(i) = plot(data.time,data.Cancer*cells(1,8)','red','LineWidth',1.2);
            else 
                p(i)=plot(data.time,data.Cancer*cells(1,8)','color',[0,0,0]+col(i),'LineWidth',1.2);
            end
            q = plot(last_follow_up,data.Cancer(pos)*cells(1,8),'b*','MarkerSize',12);
            q1 = plot(first_follow_up,data.Cancer(pos1)*cells(1,8),'m*','MarkerSize',8);
            hold off
        end
        int = yline(data{1,9}*cells(1,8)','LineStyle','--','color',[0,0.5,0],'LineWidth',1.5);
        xlabel('time (days)')
        ylabel(vars{8})
        xlim([0 T])
        leg = legend([p,int,q,q1],num2str(drug_values(1)),num2str(drug_values(2)),[num2str(drug_values(3)) '\it prescribed dose'],...
             num2str(drug_values(4)),num2str(drug_values(5)),num2str(drug_values(6)),...
             num2str(drug_values(7)),'Initial value','Last Follow-up day','First Follow-up day');
        title([tumor_status ' - cluster ' cluster])
        title(leg,[ drug_name ' doses (mg)'])
        %set(gca,'FontSize',18)
        
    end


end
