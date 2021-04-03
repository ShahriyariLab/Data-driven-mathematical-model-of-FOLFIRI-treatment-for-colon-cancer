% PlotQSPData: script for plotting the results obtained by running 
%   FOLFIRIClusterSensitivity.py
% Author: Arkadz Kirshtein, https://sites.google.com/site/akirshtein/
% (c) Shahriyari Lab https://sites.google.com/site/leilishahriyari/
%
% If using this or related code please cite 
%  Budithi, A.; Su, S.; Kirshtein,A.; Shahriyari L.
%    Data driven mathematical model of FOLFIRI treatment for colon cancer. //Cancers.
%     (Manuscript submitted for publication).

paralabels={'\lambda_{ThD5fu}', '\lambda_{5fuDTc}', '\delta_{C5fuIg}', '\delta_{5fuM}', '\alpha_{5fu}', '\delta_{5fuD}', '\delta_{5fu}', ...
                                '\delta_{CFU}', '\delta_{CFULV}', '\alpha_{IRC}', '\alpha_{IRTr}', '\delta_{CIR}', '\delta_{TrIR}', '\delta_{IR}', '\delta_{LV}', ...
                                '\alpha_{LV}'};
filenames_suffix={'_local', '_local_relative'};
plotnames={'Most sensitive parameters', 'Most relatively sensitive parameters'};

for s=1:2
	figure;
	sensvals=cell(5,2);
	for cluster=1:5
		Sensdat=dlmread(['Results/Sensitivity/V63-cluster-' num2str(cluster) '-of-5-results-sensitivity' filenames_suffix{s} '.csv']);
		for i=1:2
			subplot(5,2,2*(cluster-1)+i)
			threshold=1;
			[~,index]=sort(abs(Sensdat(1:16,i)),'descend');
			index=index(1:4);
			sensvals{cluster,i}=Sensdat(index,i);
			bar(Sensdat(index,i))
			set(gca,'fontsize',9)
			xlim([0 length(index)+1])
			xticks(1:length(index))
			xticklabels(paralabels(index))
		end
	end

	subplot(5,2,1)
	title('Sensitivity of Cancer')
	subplot(5,2,2)
	title('Sensitivity of Total Cell Count')
	for i=1:5
		subplot(5,2,2*i-1)
		ylabel(['Cluster ' num2str(i)])
		subplot(5,2,2*i)
		ylabel(['Cluster ' num2str(i)])
	end
	sgtitle(plotnames{s})
end