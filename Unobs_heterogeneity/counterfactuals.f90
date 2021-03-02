subroutine counterfactuals(params_MLE)
    use dimensions; use cadastral_maps; use simulation; use primitives
    implicit none
    double precision,dimension(par)::params_true,params_MLE
    double precision,dimension(2*P_max-1,2*P_max-1,3,3,P_max,villages)::F_true
    double precision,dimension(2*P_max-1,2,P_max,types_a,villages,unobs_types)::CCP_true
    double precision,dimension(2*P_max-1,3,P_max,types_a,unobs_types)::V_fct
    integer,dimension(plots_in_map,villages)::n_dist
    double precision,dimension(villages)::mean_N,mean_NPV,mean_budget
    integer::v_l,p_l,it
    character::end_key
    integer,parameter::nkk=21
    double precision,dimension(nkk)::pi_grid
    
    c_s_or=c_s
    c_d_or=c_d
    
    pi_grid(1)=0.0d0
    do p_l=2,nkk
        pi_grid(p_l)=pi_grid(p_l-1)+0.05d0
    end do
        
    !I want to compute the optimal tax of production giving a subsidy as a lumpsum.
    !Set a tax to production, find the lumpsum transfer that makes government transfer to be in eq.
    !Look for the optimal tax that maximizes average NPV
    
    do p_l=1,nkk;do v_l=1,1!villages
        print*,'pi exp',p_l
        print*,'village,',v_l 
        pi=pi_grid(p_l)
        c_s=c_s_or+pi*c_d_or
        c_d=c_d_or+pi*c_d_or
        if (p_l==1) then
            CCP_true(:,:,:,:,v_l,:)=0.07d0
            n_dist(:,v_l)=1
            V_fct=0.0d0
            T_g=0.0d0
        end if
1       call compute_eq_F_CCP(params_MLE,F_true(:,:,:,:,:,v_l),CCP_true(:,:,:,:,v_l,:),V_fct,n_dist(:,v_l),v_l,mean_N(v_l),mean_NPV(v_l),mean_budget(v_l))
        print*,'mean_NPV',mean_NPV(v_l)
        if (abs(mean_budget(v_l))>1.0d-3) then
            T_g=T_g+mean_budget(v_l)*0.5d0
            print*,'mean_budget(v_l)',mean_budget(v_l)
            print*,'T_g',T_g
            go to 1
        else
            print*,'eq reached for pi=',pi,'and T_g=',T_g
        end if
        if (p_l==1) then
            OPEN(UNIT=12, FILE=path_results//"counterfactuals.txt")
            write(12,'(F20.3,F20.3,F20.3,F20.3)'),pi,T_g,mean_N(v_l),mean_NPV(v_l)
            close(12)
        else
            OPEN(UNIT=12, FILE=path_results//"counterfactuals.txt",access='append')
            write(12,'(F20.3,F20.3,F20.3,F20.3)'),pi,T_g,mean_N(v_l),mean_NPV(v_l)
            close(12)
        end if        
    end do;end do
    
    !Compute NPV and mean number of wells in counterfactual

    
end subroutine
