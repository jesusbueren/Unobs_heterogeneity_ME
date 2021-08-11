
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
data_own_n=reshape(data_own_n,2,2)
fileID = fopen('counter_own_nxa.txt','r');
counter_own_n = fscanf(fileID,'%f')
fclose(fileID);
counter_own_n=reshape(counter_own_n,2,2)
fileID = fopen('modl_own_nxa.txt','r');
modl_own_n = fscanf(fileID,'%f')
fclose(fileID);
modl_own_n=reshape(modl_own_n,2,2)

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
data_uhe=reshape(data_uhe,3,2)
fileID = fopen('counter_uhe.txt','r');
counter_uhe = fscanf(fileID,'%f')
fclose(fileID);
counter_uhe=reshape(counter_uhe,3,2)
fileID = fopen('modl_uhe.txt','r');
modl_uhe = fscanf(fileID,'%f')
fclose(fileID);
modl_uhe=reshape(modl_uhe,3,2)


clrs = [0.9 0.9 0.9;0 0 0];

FS=12
set(groot,'defaultAxesTickLabelInterpreter','latex');  
for i=[2 4]
    figure(i)
    if i==2
        title("Functioning Wells in Adjacency")
        data_gr=[data_N(:) modl_N(:)] 
%         xticks([1 2 3 4 5])
%         xticklabels({'N=0','N=1','N=2','N=3','N=4'})
    elseif i==4
        title("Plots in Adjacency")
        data_gr=[data_P(2:5) modl_P(2:5)]
        xticks([1 2 3 4 ])
        xticklabels({'P=2','P=3','P=4','P=5'})
    end
%     yticks([0:0.03:0.15])
%     ylim([0 .15])
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
                errorbar(x,data_gr(:,i2), 2*sqrt((data_gr(:,i2)).*(1-data_gr(:,i2))./counter_N(:)) , 'k', 'linestyle', 'none');
            elseif i==4
                errorbar(x,data_gr(:,i2), 2*sqrt((data_gr(:,i2)).*(1-data_gr(:,i2))./counter_P(2:5)) , 'k', 'linestyle', 'none');
            end
        end
    set(hB,{'FaceColor'},{clrs(1,:),clrs(2,:)}.')
    I=legend('Data \pm 2 s.d.','Model')
    legend('boxoff')
    I.FontSize=FS
    set(gca,'FontName','Times New Roman','Fontsize',FS);
    set(gcf,'color','w')
    set(gcf,'Position',[100 100 1000 250])
    print(strcat('model_fit',num2str(i)),'-depsc')
end

set(groot,'defaultAxesTickLabelInterpreter','latex');  
i=3 
    figure(i)
    for j=1:2
        subplot(2,2,j)

            if j==1
                title("$a<1.3$",'Interpreter','latex')
            elseif j==2
                title("$1.3<a<2.3$",'Interpreter','latex')
            elseif j==3
                title("$2.3<a<4.0$",'Interpreter','latex')
            elseif j==4
                title("$a>4.0$",'Interpreter','latex')
            end
            data_gr=[data_own_n(j,:); modl_own_n(j,:)]'
%             yticks([0:0.05:0.15])
%             ylim([0 .15])
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
%             yticks([0:0.03:0.15])
%             ylim([0 .15])
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

fileID = fopen('counterfactuals_noimp.txt','r');
counterfactuals = fscanf(fileID,'%f')
fclose(fileID);

variables=5
nkk=50
villages=floor(size(counterfactuals,1)/variables/nkk) %14


counterfactuals=reshape(counterfactuals(1:variables*nkk*villages),variables,nkk,villages)

plot(mean(counterfactuals(4,:,:),3),'linewidth',2)
hold on
plot(mean(counterfactuals(5,:,:),3),'--','linewidth',2)
I=legend('Social','Private')
set(gcf,'color','w') 
        