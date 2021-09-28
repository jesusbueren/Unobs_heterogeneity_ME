
%% Counterfactuals 

close all
clear all
cd("C:\Users\jbueren\Google Drive\overdrilling\fortran\Unobs_heterogeneity_ME\Results")

fileID = fopen('counterfactuals.txt','r');
counterfactuals = fscanf(fileID,'%f')
fclose(fileID);

columns=5
lines=size(counterfactuals,1)/columns
counterfactuals=reshape(counterfactuals,columns,lines);
counterfactuals=counterfactuals';
FS=11

subplot(1,2,1)
plot(counterfactuals(1:50,1),counterfactuals(1:50,4),'linewidth',2)
xlabel('Tax per well (thousands of Rupees)')
xline(find(counterfactuals(1:50,4)== max(counterfactuals(1:50,4))),'linewidth',2)
yticks([-0.6:0.2:1])
title('Social Welfare per Acre of Land')
set(gca,'FontName','Times New Roman','Fontsize',FS);
set(gcf,'color','w')
subplot(1,2,2)
plot(counterfactuals(1:50,1),counterfactuals(1:50,3)./counterfactuals(1,3),'linewidth',2)
xlabel('Tax per well (thousands of Rupees)')
xline(find(counterfactuals(1:50,4)== max(counterfactuals(1:50,4))),'linewidth',2)
yticks([0:0.2:1])
title('Wells Relative to Benchmark')
set(gca,'FontName','Times New Roman','Fontsize',FS);
set(gcf,'color','w')
set(gcf,'Position',[100 100 750 250])
print('counterfactual','-dpdf')
print('counterfactual','-deps')


