
close all
clear all
cd("C:\Users\jbueren\Google Drive\overdrilling\fortran\Unobs_heterogeneity_ME\Results")

fileID = fopen('data_N.txt','r');
data_N = fscanf(fileID,'%f')
fclose(fileID);
fileID = fopen('counter_N.txt','r');
counter_N = fscanf(fileID,'%f')
fclose(fileID);
fileID = fopen('modl_N.txt','r');
modl_N = fscanf(fileID,'%f')
fclose(fileID);

fileID = fopen('data_own_nxa.txt','r');
data_own_n = fscanf(fileID,'%f')
fclose(fileID);
data_own_n=reshape(data_own_n,2,3)
fileID = fopen('counter_own_nxa.txt','r');
counter_own_n = fscanf(fileID,'%f')
fclose(fileID);
counter_own_n=reshape(counter_own_n,2,3)
fileID = fopen('modl_own_nxa.txt','r');
modl_own_n = fscanf(fileID,'%f')
fclose(fileID);
modl_own_n=reshape(modl_own_n,2,3)

fileID = fopen('data_P.txt','r');
data_P = fscanf(fileID,'%f')
fclose(fileID);
fileID = fopen('counter_P.txt','r');
counter_P = fscanf(fileID,'%f')
fclose(fileID);
fileID = fopen('modl_P.txt','r');
modl_P = fscanf(fileID,'%f')
fclose(fileID);

fileID = fopen('data_uhe.txt','r');
data_uhe = fscanf(fileID,'%f')
fclose(fileID);
data_uhe=reshape(data_uhe,3,3)
fileID = fopen('counter_uhe.txt','r');
counter_uhe = fscanf(fileID,'%f')
fclose(fileID);
counter_uhe=reshape(counter_uhe,3,3)
fileID = fopen('modl_uhe.txt','r');
modl_uhe = fscanf(fileID,'%f')
fclose(fileID);
modl_uhe=reshape(modl_uhe,3,3)

fileID = fopen('data_Nbar.txt','r');
data_Nbar = fscanf(fileID,'%f')
fclose(fileID);

fileID = fopen('counter_Nbar.txt','r');
counter_Nbar = fscanf(fileID,'%f')
fclose(fileID);

fileID = fopen('modl_Nbar.txt','r');
modl_Nbar = fscanf(fileID,'%f')
fclose(fileID);


fileID = fopen('data_V.txt','r');
data_V = fscanf(fileID,'%f')
fclose(fileID);
fileID = fopen('counter_V.txt','r');
counter_V = fscanf(fileID,'%f')
fclose(fileID);
fileID = fopen('modl_V.txt','r');
modl_V = fscanf(fileID,'%f')
fclose(fileID);

scatter(modl_V,data_V)
hold on
% P = polyfit(modl_V,data_V,1);
x=[0.02:0.01:0.21]
yfit = 1*x+0;
hold on;
plot(x,yfit,'r-.');
% ylim([0 0.15])
% xlim([0.05 0.1])
xlabel('Model')
ylabel('Data')

clrs = [0.9 0.9 0.9;0 0 0;0.5 0.5 0.5];

FS=9
set(groot,'defaultAxesTickLabelInterpreter','latex');  
for i=[2 4 6 7]
    figure(i)
    if i==2
        title("Functioning Wells in Adjacency")
        data_gr=[data_N(1:6) modl_N(1:6)] 
        xticks([1 2 3 4 5 6])
        xticklabels({'N=0','N=1','N=2','N=3','N=4','N=5'})
    elseif i==4
        title("Plots in Adjacency")
        data_gr=[data_P(2:6) modl_P(2:6)]
        xticks([1 2 3 4 5])
        xticklabels({'P=2','P=3','P=4','P=5','P=6'})
        
    elseif i==6
        title("Villages")
        data_gr=[data_V(:) modl_V(:)]
        xticks(1:14)
    elseif i==7
        title("N bar: ML")
        data_gr=[data_Nbar modl_Nbar]
%         xticklabels({'P=2','P=3','P=4','P=5'})
    end
    yticks([0:0.05:0.40])
    ylim([0 .4])
    hold on
    hB=bar(data_gr)
    ngroups = size(data_gr, 1);
    nbars = size(data_gr, 2);
        groupwidth = min(0.8, nbars/(nbars + 1.5));
        % Set the position of each error bar in the centre of the main bar
        % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
        for i2 = 1:1
            % Calculate center of each bar
            x = (1:ngroups) - groupwidth/2 + (2*i2-1) * groupwidth / (2*nbars);
            if i==2
                errorbar(x,data_gr(1:6,i2), 2*sqrt((data_gr(1:6,i2)).*(1-data_gr(1:6,i2))./counter_N(1:6)) , 'k', 'linestyle', 'none');
            elseif i==4
                errorbar(x,data_gr(:,i2), 2*sqrt((data_gr(:,i2)).*(1-data_gr(:,i2))./counter_P(2:6)) , 'k', 'linestyle', 'none');
            elseif i==6
                errorbar(x,data_gr(:,i2), 2*sqrt((data_gr(:,i2)).*(1-data_gr(:,i2))./counter_V(:)) , 'k', 'linestyle', 'none');
            elseif i==7
                errorbar(x,data_gr(:,i2), 2*sqrt((data_gr(:,i2)).*(1-data_gr(:,i2))./counter_Nbar) , 'k', 'linestyle', 'none');
            end
        end
    set(hB,{'FaceColor'},{clrs(1,:),clrs(2,:)}.')
    I=legend('Data \pm 2 s.d.','Model')
    legend('boxoff')
    I.FontSize=FS
    set(gca,'FontName','Times New Roman','Fontsize',FS);
    set(gcf,'color','w')
    set(gcf,'Position',[100 300 500 250])
    print(strcat('model_fit',num2str(i)),'-depsc')
end

set(groot,'defaultAxesTickLabelInterpreter','latex');  
i=3 
    figure(i)
    set(i,'position',[50    150    450    200])
    for j=1:2
        subplot(1,2,j)
            if j==1
                title("$a<1.3$",'Interpreter','latex')
            elseif j==2
                title("$1.3<a$",'Interpreter','latex')
            elseif j==3
                title("$2.3<a<4.0$",'Interpreter','latex')
            elseif j==4
                title("$a>4.0$",'Interpreter','latex')
            end
            data_gr=[data_own_n(j,:); modl_own_n(j,:)]'
            yticks([0:0.05:0.35])
            ylim([0 0.35])
            set(gca,'TickLabelInterpreter','latex')
        xticks([1 2])
        xticklabels({'$n=0$','$n=1$'})
        hold on
        hB=bar(data_gr)
        ngroups = size(data_gr, 1);
        nbars = size(data_gr, 2);
        groupwidth = min(0.8, nbars/(nbars + 1.5));
        for i2 = 1:1
            % Calculate center of each bar
            x = (1:ngroups) - groupwidth/2 + (2*i2-1) * groupwidth / (2*nbars);
            if i==3
                errorbar(x,data_gr(:,i2), 2*sqrt((data_gr(:,i2)).*(1-data_gr(:,i2))./counter_own_n(j,:)') , 'k', 'linestyle', 'none');
            elseif i==5
                errorbar(x,data_gr(:,i2), 2*sqrt((data_gr(:,i2)).*(1-data_gr(:,i2))./counter_uhe(j,:)') , 'k', 'linestyle', 'none');
            end
        end
        set(hB,{'FaceColor'},{clrs(1,:),clrs(2,:)}.')
        I=legend('Data \pm 2 s.d.','Model')
        legend('boxoff')
        I.FontSize=FS
        set(gca,'FontName','Times New Roman','Fontsize',FS);
        set(gcf,'color','w') 
%         set(gcf,'Position',[100 100 1000 500])
        print(strcat('model_fit',num2str(i)),'-depsc')
    end 


i= 5
    figure(i)
    for j=1:3
        subplot(2,2,j)
        if j==1
            title("Type I",'Interpreter','latex')
        elseif j==2
            title("Type II",'Interpreter','latex')
        elseif j==3
            title("Type III",'Interpreter','latex')
        end
            data_gr=[data_uhe(j,:); modl_uhe(j,:)]'
%             yticks([0:0.05:0.4])
%             ylim([0 .4])
             set(gca,'TickLabelInterpreter','latex')

        xticks([1 2])
        xticklabels({'$n=0$','$n=1$'})
        hold on
        hB=bar(data_gr)
        ngroups = size(data_gr, 1);
        nbars = size(data_gr, 2);
        groupwidth = min(0.8, nbars/(nbars + 1.5));
        for i2 = 1:1
            % Calculate center of each bar
            x = (1:ngroups) - groupwidth/2 + (2*i2-1) * groupwidth / (2*nbars);
            if i==3
                errorbar(x,data_gr(:,i2), 2*sqrt((data_gr(:,i2)).*(1-data_gr(:,i2))./counter_own_n(j,:)') , 'k', 'linestyle', 'none');
            elseif i==5
                errorbar(x,data_gr(:,i2), 2*sqrt((data_gr(:,i2)).*(1-data_gr(:,i2))./counter_uhe(j,:)') , 'k', 'linestyle', 'none');
            end
        end
        set(hB,{'FaceColor'},{clrs(1,:),clrs(2,:)}.')
        I=legend('Data \pm 2 s.d.','Model')
        legend('boxoff')
        I.FontSize=FS
        set(gca,'FontName','Times New Roman','Fontsize',FS);
        set(gcf,'color','w') 
        set(gcf,'Position',[100 100 1000 500])
        print(strcat('model_fit',num2str(i)),'-depsc')
    end 


%% Counterfactuals

close all
clear all
cd("C:\Users\jbueren\Google Drive\overdrilling\fortran\Unobs_heterogeneity_ME\Results")
FS=9
fileID = fopen('counterfactuals_noimp_ts.txt','r');
counterfactuals = fscanf(fileID,'%f')
fclose(fileID);

variables=5
nkk=18
villages=floor(size(counterfactuals,1)/variables/nkk) %14
beta=0.98

counterfactuals=reshape(counterfactuals(1:variables*nkk*villages),variables,nkk,villages)


 v_l=1
% A=mean(counterfactuals(4,:,:),3)
% B=mean(counterfactuals(5,:,:),3)
A=counterfactuals(4,:,v_l)
c_e=11
figure(1)
set(1,'position',[50    500    550    200])
subplot(1,2,1)
plot([0:1:nkk-1],A,'linewidth',2)
[yMax, xMax] = max(A); % xMax is an integer index 1,2,3, or 4,.....not a floating point value.
yl = ylim();
line([c_e, c_e], [yl(1), A(c_e+1)],'color',[0.5,0.5,0.5 ], 'LineWidth', 2);
xticks([0:5:nkk])
xlim([0 nkk])
set(gcf,'color','w') 
xlabel('Tax Per Well (000s Rs)')
title('Social Welfare (R/acre)')
set(gca,'FontName','Times New Roman','Fontsize',FS);
subplot(1,2,2)
A=counterfactuals(3,:,v_l)
plot([0:1:nkk-1],A,'linewidth',2)
[yMax, xMax] = max(A); % xMax is an integer index 1,2,3, or 4,.....not a floating point value.
yl = ylim();
line([c_e, c_e], [yl(1), A(c_e+1)],'color',[0.5,0.5,0.5 ], 'LineWidth', 2);
xticks([0:5:nkk])
xlim([0 nkk])
set(gcf,'color','w') 
xlabel('Tax Per Well (000s Rs)')
title('Wells per Plot')
set(gca,'FontName','Times New Roman','Fontsize',FS);
print('counterfactuals','-depsc')


%% Graphs presentation

close all
clear all
cd("C:\Users\jbueren\Google Drive\overdrilling\fortran\Unobs_heterogeneity_ME\primitives")
FS=9
TAB = readtable('flow_fail_prob_r.csv');
A=TAB{:,:}
B=reshape(A(:,4:8),3,2,18,5)
Q=[0.1;0.25;0.5;0.75;1];
figure(7)
set(7,'position',[50    150    450    200])
plot(squeeze(B(1,1,:,:))*Q*0.5+squeeze(B(1,2,:,:))*Q*0.5,'linewidth',2)
hold on
plot(squeeze(B(2,1,:,:))*Q*0.5+squeeze(B(2,2,:,:))*Q*0.5,'linewidth',2)
% plot(squeeze(B(1,1,:,:))*Q*0.5+squeeze(B(1,2,:,:))*Q*0.5,'linewidth',2)
% plot(squeeze(B(3,1,:,:))*Q*0.5+squeeze(B(3,2,:,:))*Q*0.5,'linewidth',2)
ylim([0.2 .7])
xlim([1 10])
xlabel('Number of Wells')
set(gcf,'color','w') 
set(gca,'FontName','Times New Roman','Fontsize',FS);
print('flow_pr','-depsc')

B=reshape(A(:,9),3,2,18)
figure(2)
set(2,'position',[50    150    225    200])
plot(squeeze(B(1,1,:))*0.5+squeeze(B(1,2,:))*0.5,'linewidth',2)
hold on
plot(squeeze(B(2,1,:))*0.5+squeeze(B(2,2,:))*0.5,'linewidth',2)
yticks([0:0.1:0.7])
ylim([0 0.7])
xlim([1 10])
xlabel('Number of Wells')
set(gca,'FontName','Times New Roman','Fontsize',FS);
set(gcf,'color','w') ;
print('failure_pr','-depsc')

%% Transitional Dynamics

close all
clear all
cd("C:\Users\jbueren\Google Drive\overdrilling\fortran\Unobs_heterogeneity_ME\Unobs_heterogeneity")

fileID = fopen('transitional_dynamics.txt','r');
dynamics = fscanf(fileID,'%f')
fclose(fileID);

columns=6
lines=size(dynamics,1)/columns
dynamics=reshape(dynamics,columns,lines);
dynamics=dynamics';
FS=11

T2=100
subplot(1,2,1)
plot([-5:1:T2],dynamics(16:T2+21,4),'linewidth',2)
xticks([0:20:T2])
xlim([-10 T2+5])
xlabel('Time (Years)')
yticks([-0.5:0.2:0.5])
ylim([-0.5 0.5])
title('Social Welfare per Acre of Land')
set(gca,'FontName','Times New Roman','Fontsize',FS);
set(gcf,'color','w')
subplot(1,2,2)
plot([-5:1:T2],dynamics(16:T2+21,3)./dynamics(16,3),'linewidth',2)

xlabel('Time (Years)')
yticks([0.2:0.2:1])
xticks([0:20:T2])
xlim([-10 T2+5])
ylim([0.2 1.05])
title('Wells Relative to Benchmark')
set(gca,'FontName','Times New Roman','Fontsize',FS);
set(gcf,'color','w')
set(gcf,'Position',[50    150    550    200])
cd('C:\Users\jbueren\Google Drive\overdrilling\fortran\Unobs_heterogeneity_ME\Results')
print('dynamics','-depsc')


%% Quantifying no strategic interactions

close all
clear all
cd("C:\Users\jbueren\Google Drive\overdrilling\fortran\Unobs_heterogeneity_ME\Unobs_heterogeneity")

fileID = fopen('wells_benchmark.txt','r');
benchmark = fscanf(fileID,'%f');
fclose(fileID);

fileID = fopen('wells_counterfactual.txt','r');
counterfactual = fscanf(fileID,'%f');
fclose(fileID);
EDGES=0.38:0.003:0.5

FS=11
histogram(benchmark,EDGES,'Normalization','probability')
hold on
histogram(counterfactual,EDGES,'Normalization','probability')
legend('Benchmark','No Interactions')
set(gca,'FontName','Times New Roman','Fontsize',FS);
set(gcf,'color','w')
set(gcf,'Position',[50    150    450    200])
cd('C:\Users\jbueren\Google Drive\overdrilling\fortran\Unobs_heterogeneity_ME\Results')
print('no_interactions','-depsc')


        