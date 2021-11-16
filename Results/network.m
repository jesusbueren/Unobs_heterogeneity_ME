close all
clear all
cd("C:\Users\jbueren\Google Drive\overdrilling\fortran\Unobs_heterogeneity_ME\cadastral_maps\fortran_files\")
fileID = fopen('new_map_4.txt','r');
V = fscanf(fileID,'%f');
fclose(fileID);


V_mat=reshape(V,sqrt(size(V,1)),sqrt(size(V,1)));
V_mat=V_mat-diag(diag(V_mat));
n_nodes=517
for j=1:2
    if j==1
        fileID = fopen('map_types_iid_sc_4.txt','r');
        types = fscanf(fileID,'%f')
        fclose(fileID);
    else
        fileID = fopen('map_types_99_sc_4.txt','r');
        types = fscanf(fileID,'%f')
        fclose(fileID);
    end        

node_color=zeros(size(types,1),3);
for i=1:n_nodes
    if types(i)==1
        node_color(i,:)=[1 0 0];
    elseif types(i)==2
        node_color(i,:)=[0 1 0];
    else
        node_color(i,:)=[0 0 1];
    end
end

if j==1
G=graph(V_mat(1:n_nodes,1:n_nodes))
end
figure(j)
h=plot(G,'NodeColor',node_color(1:n_nodes,:))
h.NodeLabel = {};
if j==1
title('rho=0.3')
else
title('rho=0.99')
end
end

