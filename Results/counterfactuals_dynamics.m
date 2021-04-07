
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

%% Transitional Dynamics

close all
clear all
cd("C:\Users\jbueren\Google Drive\overdrilling\fortran\Unobs_heterogeneity_ME\Results")

fileID = fopen('transitional_dynamics.txt','r');
dynamics = fscanf(fileID,'%f')
fclose(fileID);

columns=5
lines=size(dynamics,1)/columns
dynamics=reshape(dynamics,columns,lines);
dynamics=dynamics';
FS=11

subplot(1,2,1)
plot([-5:1:30],dynamics(16:51,4),'linewidth',2)
xticks([-5:5:30])
xlabel('Time (Years)')
yticks([-0.8:0.4:1])
ylim([-0.8 1])
title('Social Welfare per Acre of Land')
set(gca,'FontName','Times New Roman','Fontsize',FS);
set(gcf,'color','w')
subplot(1,2,2)
plot([-5:1:30],dynamics(16:51,3)./dynamics(16,3),'linewidth',2)
xticks([-5:5:30])
xlabel('Time (Years)')
yticks([0:0.2:1.15])
xticks([-5:5:30])
xlim([-10 35])
ylim([0.35 1.05])
title('Wells Relative to Benchmark')
set(gca,'FontName','Times New Roman','Fontsize',FS);
set(gcf,'color','w')
set(gcf,'Position',[100 100 750 250])
print('dynamics','-dpdf')
print('dynamics','-deps')