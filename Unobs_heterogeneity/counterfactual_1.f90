subroutine counterfactual_1(params_MLE)
    use dimensions; use cadastral_maps; use simulation; use primitives
    implicit none
    double precision,dimension(par)::params_true,params_MLE
    double precision,dimension(2*P_max-1,2*P_max-1,3,3,P_max,villages)::F_true
    double precision,dimension(2*P_max-1,2,P_max,types_a,villages,unobs_types)::CCP_true
    double precision,dimension(2*P_max-1,3,P_max,types_a,unobs_types)::V_fct
    integer,dimension(plots_in_map,villages)::n_dist
    double precision,dimension(villages)::mean_budget
    integer,parameter::samples=1
    double precision,dimension(samples)::mean_N,mean_NPV
    integer::v_l,p_l,it,ns
    character::end_key
    integer,parameter::nkk=8
    double precision,dimension(nkk)::fraction_grid

    
    tau=0.0d0
    
    fraction_grid(1)=1.0d0
    do p_l=2,nkk
        fraction_grid(p_l)=fraction_grid(p_l-1)-0.1d0
    end do
        
    !I want to look for the optimal NPV by selecting a given number of farmers
    do v_l=1,1!villages
        print*,'village,',v_l 
        !Run benchmark and compute NPV
        CCP_true(:,:,:,:,v_l,:)=0.07d0
        n_dist(:,v_l)=1
        V_fct=0.0d0
        T_g=0.0d0
        call compute_eq_F_CCP(params_MLE,F_true(:,:,:,:,:,v_l),CCP_true(:,:,:,:,v_l,:),V_fct,n_dist(:,v_l),v_l,mean_N(1),mean_NPV(1),mean_budget(v_l),Pr_u_X(:,:,:,:,v_l,:))
        !print*,'mean_NPV',mean_NPV(v_l)
        !OPEN(UNIT=12, FILE=path_results//"counterfactuals_1.txt")
        !    write(12,'(F20.3,F20.3,F20.3)'),0.0,mean_N(1),mean_NPV(1)
        !close(12)
        !do p_l=1,nkk
        !    pr_non_zombie(v_l)=0.01d0 !pr_non_zombie(v_l)*0.95d0
        !    do ns=1,samples
        !        print*,'p=',p_l,' ns=',ns
        !        call load_cadastral_maps()
        !        call compute_eq_F_CCP(params_MLE,F_true(:,:,:,:,:,v_l),CCP_true(:,:,:,:,v_l,:),V_fct,n_dist(:,v_l),v_l,mean_N(ns),mean_NPV(ns),mean_budget(v_l),Pr_u_X(:,:,:,:,v_l,:))
        !        print*,'mean_NPV: ',mean_NPV(ns)
        !    end do
        !    OPEN(UNIT=12, FILE=path_results//"counterfactuals_1.txt",access='append')
        !    write(12,'(F20.3,F20.3,F20.3)'),fraction_grid(p_l),mean_N(maxloc(mean_NPV)),maxval(mean_NPV)
        !    close(12)
        !end do
       
    end do
    
    !Compute NPV and mean number of wells in counterfactual

    
end subroutine
    


    