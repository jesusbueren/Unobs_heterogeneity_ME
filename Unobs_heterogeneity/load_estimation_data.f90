subroutine load_estimation_data
    use simulation
    implicit none
    double precision,dimension(21,T_sim,plots_i)::data_csv
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
    UHE_type=0.25d0!data_csv(16:18,1,:)
    drilling_it(:,:,1)=data_csv(19,:,:) 
    impute_i=data_csv(20,1,:)
    can_be_zombie_i=data_csv(21,1,:)
    !print*,'no zombie in estimation data'
    !impute_i=can_be_zombie_i
    can_be_zombie_i=0
    

    !UHE_type(selected_type,:)=1.0d0
    
    !Modal value of unobserved heterogeneity
    do i_l=1,plots_i
        modal_UHE_type(i_l)=maxloc(UHE_type(:,i_l),1)
    end do
    
    !Modal number of functionning wells
    do i_l=1,plots_i;do t_l=1,T_sim
        modal_N(t_l,i_l)=maxloc(Pr_N_data(:,t_l,i_l),1)-(n_data(t_l,i_l)-1)
    end do;end do
    
    
    
end subroutine
