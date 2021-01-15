module dimensions
    implicit none
    integer,parameter::P_max=6 ! Set the maximum number of plots in an adjacency
    integer,parameter::K=5,par=4,M=2,types_a=4 !K: points of support of flow; M:types of moonzoons; type_a: types of areas
end  
    
module cadastral_maps
    use dimensions
    implicit none
    integer,parameter::plots_in_map=1909,villages=14,unobs_types=4
    integer,parameter,dimension(villages)::plots_v=(/1794,302,912,517,292,535,939,637,405,837,973,1844,443,1909/) !plots in each village
    integer,dimension(plots_in_map,plots_in_map,villages)::neighbors_map
    integer,dimension(plots_in_map,2,villages)::PA_type !number of neighbors, a type for each plot in the map
    integer,dimension(plots_in_map,P_max,villages)::neighbors !identify neighbors for each plot in the map
    character(len=103)::file_map="C:\Users\jbueren\Google Drive\overdrilling\fortran\Unobs_heterogeneity_ME\cadastral_maps\fortran_files\"
end
    
module primitives
    use dimensions; use cadastral_maps
    implicit none
    character(len=85)::path_primitives="C:\Users\jbueren\Google Drive\overdrilling\fortran\Unobs_heterogeneity_ME\primitives\"
    !q: Discharge distribution (flow for each point of the support)
    double precision,dimension(K,1)::q
    !PI_k: Discharge pr. (first position indicates 1 well) depending on n0 of wells
    double precision,dimension(2*P_max,K,M,unobs_types)::PI_k
    !PI_m: unconditional pr of high and low monzoon
    double precision,dimension(M,villages)::PI_m
    !PI_f: Pr. of failure (first position indicates 1 well; last position indicates all plots have 2 wells)
    double precision,dimension(2*P_max,M,unobs_types)::PI_fm
    double precision,dimension(2*P_max-1,3,P_max,villages,unobs_types)::PI_f_v
    !PI_s: Pr. of success (first position indicates no wells; last position indicates all plots 2 wells but one with one well)
    double precision,dimension(2*P_max,villages)::PI_s
    double precision,dimension(2*P_max-1,3,P_max,villages)::PI_s_v
    !c_d: fixed cost of failing to drill;c_s: fixed cost of succeeding to drill
    double precision,parameter::c_s=66.4d0,beta=0.8d0,c_d=27.2d0
    !extreme value distribution shocks
    double precision,parameter::gamma=0.577215664901533d0
    double precision::rho
    !area of plots
    double precision,dimension(types_a)::area=(/1.0d0,2.0d0,3.0d0,5.1d0/)!(/0.394d0,1.279d0,3.157d0,8.706d0/)
    !pr of unobserved heterogeneity type
    double precision,dimension(unobs_types)::pr_unobs_t=(/0.216d0*0.168d0,0.216d0*0.832d0,0.784d0*0.168d0,0.784d0*0.832d0/)
    
end
    
module simulation
use cadastral_maps
    implicit none
    character(len=79)::path_estimation="C:\Users\jbueren\Google Drive\overdrilling\fortran\Unobs_heterogeneity_ME\data\"
    character(len=82)::path_results="C:\Users\jbueren\Google Drive\overdrilling\fortran\Unobs_heterogeneity_ME\Results\"
    !maximum number of functioning wells seen in the data drill_export_.xls
    integer,parameter::max_NFW=10,simulations=1
    !Parameters of the simulated panel-> I: number of indv; T: number of periods
    integer,parameter::T_sim=5,plots_i=1052
    ! dec_it: drilling decision
    double precision,dimension(2*P_max-1,2*P_max-1,3,3,P_max,villages)::F_est
    double precision,dimension(2*P_max-1,2,P_max,types_a,villages,unobs_types)::CCP_est
    
    !Data
    integer,dimension(plots_i)::V_type,P_type,A_type
    double precision,dimension(unobs_types,plots_i)::UHE_type
    integer,dimension(plots_i)::modal_UHE_type
    integer,dimension(T_sim,plots_i,simulations)::drilling_it
    integer,dimension(T_sim,plots_i)::n_data !number of wells in reference plot
    double precision,dimension(max_NFW+1,T_sim,plots_i)::Pr_N_data !pr of number of functioning wells in the adjacency
    integer,dimension(T_sim,plots_i)::modal_N !pr of number of functioning wells in the adjacency
    
end module    
    