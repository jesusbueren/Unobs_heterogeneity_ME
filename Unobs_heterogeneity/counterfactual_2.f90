subroutine counterfactual_2(params_MLE)
    use dimensions; use cadastral_maps; use simulation; use primitives
    implicit none
    double precision,dimension(par)::params_true,params_MLE
    double precision,dimension(2*P_max-1,2*P_max-1,3,3,P_max,villages,unobs_types)::F_true
    double precision,dimension(2*P_max-1,2,P_max,types_a,villages,unobs_types)::CCP_true
    double precision,dimension(2*P_max-1,3,P_max,types_a,unobs_types)::V_fct,V_social
    integer,dimension(plots_in_map,villages)::n_dist
    double precision,dimension(villages)::mean_N,social_output,private_output
    integer::v_l,p_l,it
    character::end_key
    integer,parameter::nkk=20
    double precision,dimension(nkk)::tau_grid
    
    rho=params_MLE(3)
    print*,'p2',params_MLE
    tau_grid(1)=0.0d0
    do p_l=2,nkk
        tau_grid(p_l)=tau_grid(p_l-1)+1.0d0
    end do
    
        
    !I want to compute the optimal tax of production giving a subsidy as a lumpsum.
    !Set a tax to production, find the lumpsum transfer that makes government transfer to be in eq.
    !Look for the optimal tax that maximizes average NPV
    tau=0.0d0
    do v_l=1,1;do p_l=1,nkk
        print*,'exp',p_l
        print*,'village,',v_l 
        tau=tau_grid(p_l)
        if (p_l==1) then
            CCP_true(:,:,:,:,v_l,:)=0.15d0
            n_dist(:,v_l)=1
            V_fct=0.0d0
            V_social=0.0d0
            !T_g=0.0d0
        end if
!        !If needed to compute the equilibrium tax
!1       call compute_eq_F_CCP(params_MLE,F_true(:,:,:,:,:,v_l),CCP_true(:,:,:,:,v_l,:),V_fct,n_dist(:,v_l),v_l,mean_N(v_l),mean_NPV(v_l),mean_budget(v_l))
!        print*,'mean_NPV',mean_NPV(v_l)
!        if (abs(mean_budget(v_l))>1.0d-2) then
!            T_g=T_g+mean_budget(v_l)*0.5d0
!            print*,'mean_budget(v_l)',mean_budget(v_l)
!            print*,'T_g',T_g
!            go to 1
!        else
!            print*,'eq reached for tau=',tau,'and T_g=',T_g
!        end if
        call compute_eq_F_CCP(params_MLE,F_true(:,:,:,:,:,v_l,:),CCP_true(:,:,:,:,v_l,:),V_fct,V_social,n_dist(:,v_l),v_l,mean_N(v_l),social_output(v_l),private_output(v_l),Pr_u_X(:,:,:,:,v_l,:))
        if (p_l==1 .and. v_l==1) then
            OPEN(UNIT=12, FILE=path_results//"counterfactuals_noimp_ts.txt")
            write(12,'(F10.3,I4,F10.3,F10.3,F10.3)'),tau,v_l,mean_N(v_l),social_output(v_l),private_output(v_l)
            close(12)
        else
            OPEN(UNIT=12, FILE=path_results//"counterfactuals_noimp_ts.txt",access='append')
            write(12,'(F10.3,I4,F10.3,F10.3,F10.3)'),tau,v_l,mean_N(v_l),social_output(v_l),private_output(v_l)
            close(12)
        end if        
    end do;end do


    
end subroutine
