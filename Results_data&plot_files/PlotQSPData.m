%Plot dynamics - Figure 3

% reading parsed data
load('Results\All_Dynamic\V63-results-all-data.mat')

% total time = start-treatment + treatment + observation time 
T=1000;

figure;
cmap=lines(5);
for i=1:9
    subplot(3,6,i);
    hold on
    for cluster=1:5
        plot(data{cluster}(:,1),data{cluster}(:,i+1)*cells(cluster,i),'color',cmap(cluster,:),'LineWidth',1);
    end
    xlabel('time (days)','FontSize',15)
    ylabel(vars{i},'FontSize',15)
    xlim([0 T])
    hold off
end
%drugs
for i=15:17
    subplot(3,6,i-5);
    hold on
    for cluster=1:5
        plot(data{cluster}(:,1),data{cluster}(:,i+1),'color',cmap(cluster,:),'LineWidth',1);
    end
    xlabel('time (days)','FontSize',15)
    ylabel(vars{i},'FontSize',15)
    xlim([0 T])
    hold off
end
subplot(3,6,18);
hold on
for cluster=1:5
    plot(data{cluster}(:,1),data{cluster}(:,2:10)*cells(cluster,lmod)','color',cmap(cluster,:),'LineWidth',1);
end
xlabel('time (days)','FontSize',15)
ylabel(vars{end},'FontSize',15)
xlim([0 T])
legend({'cluster 1','cluster 2','cluster 3','cluster 4','cluster 5'},'LineWidth',1,'FontSize',12)
hold off


%cytokines
for i=10:14
    subplot(3,6,i+3);
    hold on
    for cluster=1:5
        plot(data{cluster}(:,1),data{cluster}(:,i+1)*cells(cluster,i),'color',cmap(cluster,:),'LineWidth',1);
    end
    xlabel('time (days)','FontSize',15)
    ylabel(vars{i},'FontSize',15)
    xlim([0 T])
    hold off
end
