
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
fileID = fopen('modl_own_nxa.txt','r');
modl_own_n = fscanf(fileID,'%f')
fclose(fileID);

fileID = fopen('data_P.txt','r');
data_P = fscanf(fileID,'%f')
fclose(fileID);
fileID = fopen('modl_P.txt','r');
modl_P = fscanf(fileID,'%f')
fclose(fileID);

fileID = fopen('data_uhe.txt','r');
data_uhe = fscanf(fileID,'%f')
fclose(fileID);

fileID = fopen('modl_uhe.txt','r');
modl_uhe = fscanf(fileID,'%f')
fclose(fileID);



clrs = [0 0 0; 0.9 0.9 0.9];

FS=12
set(groot,'defaultAxesTickLabelInterpreter','latex');  
for i=2:5
    figure(i)
    if i==2
        title("Functioning Wells in Adjacency")
        data_gr=[data_N(1:5) modl_N(1:5)]
        xticks([1 2 3 4 5])
        xticklabels({'N=0','N=1','N=2','N=3','N=4'})
    elseif i==3
        title("Functioning Wells in Reference Plot x Area")
        data_gr=[data_own_n modl_own_n]
        xticks([1:8])
        xticklabels({'$$\begin{array}{c}a<1.3 \\ n=0 \end{array}$$',...
                     '$$\begin{array}{c}1.3<a<2.3 \\ n=0 \end{array}$$', ...
                     '$$\begin{array}{c}2.3<a<4.0 \\ n=0 \end{array}$$', ...
                     '$$\begin{array}{c}a>4.0 \\ n=0 \end{array}$$', ...
                     '$$\begin{array}{c}a<1.3 \\ n=1 \end{array}$$',...
                     '$$\begin{array}{c}1.3<a<2.3 \\ n=1 \end{array}$$', ...
                     '$$\begin{array}{c}2.3<a<4.0 \\ n=1 \end{array}$$', ...
                     '$$\begin{array}{c}a>4.0 \\ n=1 \end{array}$$'})
         set(gca,'TickLabelInterpreter','latex')
    elseif i==4
        title("Plots in Adjacency")
        data_gr=[data_P(3:6) modl_P(3:6)]
        xticks([1 2 3 4])
        xticklabels({'P=3','P=4','P=5','P=6'})
    elseif i==5
        title("Unobserved Heterogeneity")
        data_gr=[data_uhe modl_uhe]
        xticks([1:8])
        xticklabels({'$$\begin{array}{c}Low Fl, High Fa \\ n=0 \end{array}$$',...
                     '$$\begin{array}{c}Low Fl, Low Fa \\ n=0 \end{array}$$', ...
                     '$$\begin{array}{c}High Fl, High Fa \\ n=0 \end{array}$$', ...
                     '$$\begin{array}{c}High Fl, Low Fa  \\ n=0 \end{array}$$', ...
                     '$$\begin{array}{c}Low Fl, High Fa \\ n=1 \end{array}$$',...
                     '$$\begin{array}{c}Low Fl, Low Fa  \\ n=1 \end{array}$$', ...
                     '$$\begin{array}{c}High Fl, High Fa  \\ n=1 \end{array}$$', ...
                     '$$\begin{array}{c}High Fl, Low Fa  \\ n=1 \end{array}$$'})
         set(gca,'TickLabelInterpreter','latex')
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
    set(gcf,'Position',[100 100 1000 500])
    print(strcat('model_fit',num2str(i)),'-depsc')
end
    
        