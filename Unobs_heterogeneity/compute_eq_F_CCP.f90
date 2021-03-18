subroutine compute_eq_F_CCP(params,F,CCP_mid,V_fct,n_initial,v_l,mean_N,social_output,private_output)
    use cadastral_maps; use dimensions; use primitives
    implicit none
    double precision,dimension(par),intent(in)::params
    double precision,dimension(2*P_max-1,2*P_max-1,3,3,P_max),intent(out)::F
    double precision,dimension(2*P_max-1,2,P_max,types_a,unobs_types),intent(inout)::CCP_mid
    integer,dimension(plots_in_map,1),intent(inout)::n_initial
    double precision,dimension(2*P_max-1,3,P_max,types_a,unobs_types),intent(inout)::V_fct
    integer,intent(in)::v_l    
    double precision,intent(out)::mean_N,social_output,private_output
    double precision,dimension(2*P_max-1,2,P_max,types_a,unobs_types)::CCP_old,CCP
    
    double precision,dimension(2*P_max-1,3,P_max,types_a,villages,unobs_types)::Ef_v !Ef_v: expected productivity
    double precision::dist
    integer::p_l,a_l,n_l,P_l2,ind,counter_all,counter_bad,u_l
    integer(8),dimension(2*P_max-1,3,3,P_max)::iterations
    character::pause_k
    
    
    !Set scale parameter Gumbel distribution of shocks
    rho=params(4)

    !Compute expected productivity 
    do u_l=1,unobs_types;do a_l=1,types_a
        call expected_productivity(params(1:3),area(a_l),Ef_v(:,:,:,a_l,v_l,u_l),v_l,u_l)
    end do;end do

    !Generate beliefs consitent with CCP
    F=1.0d0
    CCP=CCP_mid
!   print*,'generating beliefs'
1    n_initial=1
    call generate_beliefs(CCP_mid,V_fct,Ef_v(:,:,:,:,v_l,:),n_initial,F,v_l,iterations,mean_N,social_output,private_output)
    
    !For each plot type obtain a new CCP given beliefs
    !print*,'policy step'
    CCP_old=CCP
    dist=0.0d0
    counter_bad=0
    counter_all=0
    do P_l=1,P_max; do a_l=1,types_a; do u_l=1,unobs_types
        call policy_fct_it(Ef_v(1:2*P_l-1,:,P_l,a_l,v_l,u_l)&
                            ,F(1:2*P_l-1,1:2*P_l-1,:,:,P_l) &
                            ,P_l &
                            ,CCP_old(1:2*P_l-1,:,P_l,a_l,u_l),CCP(1:2*P_l-1,:,P_l,a_l,u_l),v_l,u_l &
                            ,V_fct(1:2*P_l-1,:,P_l,a_l,u_l))
        !print*,''
        !print*,'V_fct',V_fct(1:2*P_l-1,:,P_l,a_l,u_l)
        !print*,'CCP',CCP(1:2*P_l-1,:,P_l,a_l,u_l)
        !call value_fct_it(Ef_v(1:2*P_l-1,:,P_l,a_l,v_l,u_l)&
        !                    ,F(1:2*P_l-1,1:2*P_l-1,:,:,P_l) &
        !                    ,P_l &
        !                    ,CCP(1:2*P_l-1,:,P_l,a_l,u_l),v_l,u_l &
        !                    ,V_fct(1:2*P_l-1,:,P_l,a_l,u_l))
        !print*,''
        !print*,'V_fct',V_fct(1:2*P_l-1,:,P_l,a_l,u_l)
        !print*,'CCP',CCP(1:2*P_l-1,:,P_l,a_l,u_l)
    end do; end do;end do
    !V_fct(1,3,2:7,1,1)
   
    !V_fct(1,3,1,:,4)
    !P_l2=P_max
    !do u_l=1,unobs_types;do a_l=1,types_a; do n_l=1,1;do ind=1,2*P_l2-1; 
    !    if (abs(CCP_old(ind,n_l,P_l2,a_l,u_l)-CCP(ind,n_l,P_l2,a_l,u_l))>0.02) then
    !        print*,'U_T',u_l,'N',n_l-2+ind,' ;A',a_l,'n_l',n_l,'ind',ind
    !        print '(F8.4,F8.4,I12,I12)',real(CCP(ind,n_l,P_l2,a_l,u_l)),real(CCP_old(ind,n_l,P_l2,a_l,u_l)),iterations(ind,n_l,1,P_l2),iterations(ind,n_l,2,P_l2)
    !        print*,''
    !    end if
    !end do; end do;end do;end do
    
    P_l2=P_max
    dist=0.0
    do P_l=2,P_max; do n_l=1,2;do ind=1,2*P_l-1; 
        dist=dist+dble(sum(iterations(ind,n_l,1:3,P_l)))/dble(sum(iterations(:,1:2,1:3,:)))*sum(abs(CCP_old(ind,n_l,P_l,:,:)-CCP(ind,n_l,P_l,:,:)))/dble(types_a)/dble(unobs_types)
    end do;end do; end do
    !print*,'village',v_l
    print*,'dist CCP',dist,'social_output',social_output
    
    !New guess of the ccp is half way through
    CCP_mid=CCP*0.5d0+CCP_old*0.5d0
    
    !print*,'npv',mean_NPV,'mean wells',mean_N
    
    !print*,'press any key to continue'
    !read*,pause_k
    if (dist>1.0d-4 ) then !dist>1d-3 
        go to 1 
    end if
    
    call generate_beliefs(CCP_mid,V_fct,Ef_v(:,:,:,:,v_l,:),n_initial,F,v_l,iterations,mean_N,social_output,private_output)
    !print*,'dist CCP',dist,'social_output',social_output
end subroutine