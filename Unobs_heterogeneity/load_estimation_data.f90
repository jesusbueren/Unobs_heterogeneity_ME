subroutine load_estimation_data
    use simulation
    implicit none
    double precision,dimension(20,T_sim,plots_i)::data_csv
    integer::i_l,t_l
    double precision,dimension(types_a)::moment_a
    double precision,dimension(max_NFW+1)::moment_N
    double precision,dimension(2)::moment_own_n
    double precision,dimension(unobs_types)::moment_uhe
    double precision,dimension(P_max)::moment_P
        
    
    OPEN(UNIT=12, FILE=path_estimation//"drill_export_r.csv")
    read(12,*),data_csv
    close(12)
    
    V_type=data_csv(1,1,:)
    P_type=data_csv(2,1,:)
    A_type=data_csv(3,1,:)
    n_data=data_csv(4,:,:)+1
    Pr_N_data=data_csv(5:15,:,:)
    UHE_type=data_csv(16:19,1,:)
    drilling_it(:,:,1)=data_csv(20,:,:)    
    
    !Modal value of unobserved heterogeneity
    do i_l=1,plots_i
        modal_UHE_type(i_l)=maxloc(UHE_type(:,i_l),1)
    end do
    
    !Modal number of functionning wells
    do i_l=1,plots_i;do t_l=1,T_sim
        modal_N(t_l,i_l)=maxloc(Pr_N_data(:,t_l,i_l),1)
    end do;end do
    
    
    
end subroutine
