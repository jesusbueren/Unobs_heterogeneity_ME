close all
clear all
cd("C:\Users\jbueren\Google Drive\overdrilling\fortran\Unobs_heterogeneity_ME\Results")

fileID = fopen('bootstrapped_parameters.txt','r');
se = fscanf(fileID,'%f')
fclose(fileID);

se=reshape(se,4,size(se,1)/4)

for i=1:3
    subplot(2,2,i)
    histogram(se(i,:),'Normalization','probability')
    if i==1
        title(['scale=',num2str(round(se(i,1),2)),' (',num2str(round(std(se(i,:)),2)),')'])
    elseif i==2
        title(['flow=',num2str(round(se(i,1),2)),' (',num2str(round(std(se(i,:)),2)),')'])
    elseif i==3
        title(['error var=',num2str(round(se(i,1),2)),' (',num2str(round(std(se(i,:)),2)),')'])
    end
end
% 
% figure(2)
% histogram(se(4,:),'Normalization','probability')