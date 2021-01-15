
close all
clear all
cd("C:\Users\jbueren\Google Drive\overdrilling\fortran\Unobs_heterogeneity_ME\Results")

fileID = fopen('data_N.txt','r');
data_N = fscanf(fileID,'%f')
fclose(fileID);
fileID = fopen('modl_N.txt','r');
modl_N = fscanf(fileID,'%f')
fclose(fileID);

fileID = fopen('data_own_nxa.txt','r');
data_own_n = fscanf(fileID,'%f')
fclose(fileID);
data_own_n=reshape(data_own_n,4,2)
fileID = fopen('modl_own_nxa.txt','r');
modl_own_n = fscanf(fileID,'%f')
fclose(fileID);
modl_own_n=reshape(modl_own_n,4,2)


fileID = fopen('data_P.txt','r');
data_P = fscanf(fileID,'%f')
fclose(fileID);
fileID = fopen('modl_P.txt','r');
modl_P = fscanf(fileID,'%f')
fclose(fileID);

fileID = fopen('data_uhe.txt','r');
data_uhe = fscanf(fileID,'%f')
fclose(fileID);
data_uhe=reshape(data_uhe,4,2)
fileID = fopen('modl_uhe.txt','r');
modl_uhe = fscanf(fileID,'%f')
fclose(fileID);
modl_uhe=reshape(modl_uhe,4,2)


clrs = [0 0 0; 0.9 0.9 0.9];

FS=12
set(groot,'defaultAxesTickLabelInterpreter','latex');  
for i=[2 4]
    figure(i)
    if i==2
        title("Functioning Wells in Adjacency")
        data_gr=[data_N(1:5) modl_N(1:5)]
        xticks([1 2 3 4 5])
        xticklabels({'N=0','N=1','N=2','N=3','N=4'})
    elseif i==4
        title("Plots in Adjacency")
        data_gr=[data_P(3:6) modl_P(3:6)]
        xticks([1 2 3 4])
        xticklabels({'P=3','P=4','P=5','P=6'})
    end
    hold on
    hB=bar(data_gr)
    set(hB,{'FaceColor'},{clrs(1,:),clrs(2,:)}.')
    I=legend('Data','Model')
    legend('boxoff')
    I.FontSize=FS
    set(gca,'FontName','Times New Roman','Fontsize',FS);
    set(gcf,'color','w')
%     ylim([0 0.1])
    if i==5
%         ylim([0 0.2])
    end 
    set(gcf,'Position',[100 100 700 250])
    print(strcat('model_fit',num2str(i)),'-depsc')
end


FS=12
set(groot,'defaultAxesTickLabelInterpreter','latex');  
for i=[3 5]
    figure(i)
    for j=1:4
    subplot(2,2,j)
    if i==3
        if j==1
            title("$a<1.3$", 'interpreter', 'latex')
        elseif j==2
            title("$1.3<a<2.3$", 'interpreter', 'latex')
        elseif j==3
            title("$2.3<a<4.0$", 'interpreter', 'latex')
        elseif j==4
            title("$a>4.0$", 'interpreter', 'latex')
        end
        data_gr=[data_own_n(j,:); modl_own_n(j,:)]'
        xticks([1:8])
        xticklabels({'n=0',...
                     'n=1'})
         set(gca,'TickLabelInterpreter','latex')
         ylim([0 0.15])
    elseif i==5
        if j==1
            title("Low Flow, High Failure")
        elseif j==2
            title("Low Flow, Low Failure")
        elseif j==3
            title("High Flow, High Failure")
        elseif j==4
            title("High Flow, Low Failure")
        end
        data_gr=[data_uhe(j,:); modl_uhe(j,:)]'
        xticks([1:8])
        xticklabels({'n=0',...
                     'n=1'})
         set(gca,'TickLabelInterpreter','latex')
         ylim([0 0.3])
    end
    hold on
    hB=bar(data_gr)
    set(hB,{'FaceColor'},{clrs(1,:),clrs(2,:)}.')
    I=legend('Data','Model')
    legend('boxoff')
    I.FontSize=FS
    set(gca,'FontName','Times New Roman','Fontsize',FS);
    set(gcf,'color','w')
    if i==5
%         ylim([0 0.2])
    end 
    set(gcf,'Position',[100 100 700 500])
    print(strcat('model_fit',num2str(i)),'-depsc')
    end
end
    
        