subroutine generate_beliefs(CCP,n_initial,F_new,v_l,iterations)
    use cadastral_maps; use primitives
    implicit none
    double precision,dimension(2*P_max-1,2,P_max,types_a,unobs_types),intent(in)::CCP
    integer,dimension(plots_in_map,1),intent(inout)::n_initial
    double precision,dimension(2*P_max-1,2*P_max-1,3,3,P_max),intent(out)::F_new
    integer,intent(in)::v_l
    integer(8),dimension(2*P_max-1,3,3,P_max),intent(out)::iterations
    integer,parameter::T=1000!00
    integer,dimension(plots_in_map,3)::state,state_old
    integer::i_l,j_l,t_l,ind,N_all,n_l,P,A,P_l,n_l2,it,m_l,it_min
    double precision::u_d,u_s,u_f,u_m
    integer(8),dimension(2*P_max-1,2*P_max-1,3,3,P_max)::beliefs_c
    integer,dimension(1)::seed=123
    integer,parameter::burn_t=100
    double precision,dimension(2*P_max-1,2*P_max-1,3,3,P_max)::F
    integer,dimension(plots_v(v_l))::unobs_types_i
    double precision,dimension(P_max)::dist
    character::continue_k
    
    !Call seed number
    call random_seed(PUT=seed)
    
    !Generate permanent unobserved heterogeneity type
    do i_l=1,plots_v(v_l)
        call RANDOM_NUMBER(u_m)
        if (u_m<pr_unobs_t(1)) then
            unobs_types_i(i_l)=1
        elseif (u_m<sum(pr_unobs_t(1:2))) then
            unobs_types_i(i_l)=2
        elseif (u_m<sum(pr_unobs_t(1:3))) then
            unobs_types_i(i_l)=3
        else
            unobs_types_i(i_l)=4
        end if
    end do
    
    beliefs_c=0
    iterations=0
    F_new=-9.0d0
    it=0
    
    !Store the state for each plot and simulate decision to drill
    !OPEN(UNIT=12, FILE="F_mean1.txt")
    !OPEN(UNIT=13, FILE="F_mean2.txt")
    do t_l=1,T-1;   
        !print*,t_l
        !simulate monsoon next period
        call RANDOM_NUMBER(u_m)
        if (u_m<PI_m(1,v_l))then
            m_l=1
        else
            m_l=2
        end if
        beliefs_c=0
        state(:,1)=n_initial(:,1)
        !print*,'t_l',t_l,'av number of wells per plot',real(sum(n_initial(1:plots_v(v_l),1))-plots_v(v_l))/real(plots_v(v_l))
        do i_l=1,plots_v(v_l)
            N_all=1 !Indicates the number of wells in the adjacency
            !Loop over all neighbors
            do j_l=1,PA_type(i_l,1,v_l) !PA_type(i_l,1) stores the number of plots in the adjacency
                if (state(neighbors(i_l,j_l,v_l),1)==2)  then
                    N_all=N_all+1 !number of wells (there is one well)
                elseif (state(neighbors(i_l,j_l,v_l),1)==3)  then
                    N_all=N_all+2 !number of wells (there is two wells)
                end if
            end do
            state(i_l,2)=N_all !second column in state: number of plots with one well
            n_l=state(i_l,1) !number of well in reference plot

            P=PA_type(i_l,1,v_l) !number of plots in the adjacency
            A=PA_type(i_l,2,v_l) !area of the reference plot
            
            !Locate position in the state space wrt to the CCP, PI_s_v and ,PI_f_v
            if (n_l==1) then
                ind=N_all 
            elseif (n_l==2) then
                ind=N_all-1
            elseif (n_l==3) then
                ind=N_all-2
            else
                print*,'error generating beliefs'
            end if         
            state(i_l,3)=ind
            !Count transitions (in the first iteration state_old is undefined: no problem)
            if (t_l>=burn_t) then
                beliefs_c(state_old(i_l,3),state(i_l,3),state_old(i_l,1),state(i_l,1),P)=&
                beliefs_c(state_old(i_l,3),state(i_l,3),state_old(i_l,1),state(i_l,1),P)+1
            end if
            !Well drilling decision and failures/successes
            if (n_l==1) then !no well
                call RANDOM_NUMBER(u_d)
                if (u_d<CCP(ind,n_l,P,A,unobs_types_i(i_l))) then !decides to drill
                    call RANDOM_NUMBER(u_s)
                    if (u_s<PI_s_v(ind,n_l,P,v_l)) then !successful attempt
                        n_initial(i_l,1)=n_l+1
                    else !unsuccessful attempt
                        n_initial(i_l,1)=n_l
                    end if
                else !decides not to drill
                    n_initial(i_l,1)=n_l
                end if
            elseif (n_l==2) then !one well
                call RANDOM_NUMBER(u_d)
                if (u_d<CCP(ind,n_l,P,A,unobs_types_i(i_l))) then !decides to drill
                    call RANDOM_NUMBER(u_s)
                    if (u_s<PI_s_v(ind,n_l,P,v_l)) then !successful attempt
                        call RANDOM_NUMBER(u_f)
                        if (u_f<PI_fm(N_all-1,m_l,unobs_types_i(i_l))) then !failure of the previous well
                            n_initial(i_l,1)=n_l
                        else
                            n_initial(i_l,1)=n_l+1
                        end if
                    else !unsuccessful attempt
                        call RANDOM_NUMBER(u_f)
                        if (u_f<PI_fm(N_all-1,m_l,unobs_types_i(i_l))) then !failure of the previous well
                            n_initial(i_l,1)=n_l-1
                        else
                            n_initial(i_l,1)=n_l
                        end if
                    end if
                else !decides not to drill
                    call RANDOM_NUMBER(u_f)
                    if (u_f<PI_fm(N_all-1,m_l,unobs_types_i(i_l))) then !failure of the previous well
                        n_initial(i_l,1)=n_l-1
                    else
                        n_initial(i_l,1)=n_l
                    end if 
                end if 
            elseif(n_l==3) then !two wells
                call RANDOM_NUMBER(u_f)
                if (u_f<PI_fm(N_all-1,m_l,unobs_types_i(i_l))**2) then !failure of the two wells
                    n_initial(i_l,1)=n_l-2
                elseif (u_f>(1.0d0-PI_fm(N_all-1,m_l,unobs_types_i(i_l)))**2) then !failure of none
                    n_initial(i_l,1)=n_l
                else !failure of one
                    n_initial(i_l,1)=n_l-1
                end if 
            else
                print*,'error in gen beliefs 2'
            end if
        end do
        !Store current state
        state_old=state
        !Compute beliefs
        if (t_l>burn_t+1) then
            F=0.0d0
            do P_l=2,P_max; do ind=1,2*P_l-1; do n_l=1,3; do n_l2=1,3
                    if (n_l==1 .and. n_l2==3) then
                        F(ind,:,n_l,n_l2,P_l)=-9.0d0
                    elseif ((sum(beliefs_c(ind,:,n_l,n_l2,P_l)))==0) then
                        F(ind,1:2*P_l-1,n_l,n_l2,P_l)=-9.0d0
                    else
                        F(ind,:,n_l,n_l2,P_l)=dble(beliefs_c(ind,:,n_l,n_l2,P_l))/dble(sum(beliefs_c(ind,:,n_l,n_l2,P_l)))
                        iterations(ind,n_l,n_l2,P_l)=iterations(ind,n_l,n_l2,P_l)+sum(beliefs_c(ind,:,n_l,n_l2,P_l))
                        F_new(ind,:,n_l,n_l2,P_l)=dble(sum(beliefs_c(ind,:,n_l,n_l2,P_l)))/dble(iterations(ind,n_l,n_l2,P_l))*F(ind,:,n_l,n_l2,P_l) +&
                                                    dble(iterations(ind,n_l,n_l2,P_l)-sum(beliefs_c(ind,:,n_l,n_l2,P_l)))/dble(iterations(ind,n_l,n_l2,P_l))*F_new(ind,:,n_l,n_l2,P_l) 
                        if (minval(F_new(ind,:,n_l,n_l2,P_l))<0)then
                            print*,''
                        end if
                    end if
            end do;end do;end do;end do
            !if (t_l>=burn_t+2) then
            !    write(12,*),F_new(1,:,1,2,5)
            !    write(13,*),F_new(2,:,1,2,5)
            !end if
            it=0
        end if
        it=it+1
    end do
    !close(12)
    !close(13)

!F_new(2,1:3,3,1,2) iterations(2,3,1,2)
    ! In case I don't have observations for a given state, I consider that the transition pr 
    ! is the same for all possible future states
    do P_l=2,P_max; do ind=1,2*P_l-1; do n_l=1,3; do n_l2=1,min(n_l+1,3)
        it_min=300000
        if (iterations(ind,n_l,n_l2,P_l)==0) then
            F_new(ind,1:2*P_l-1,n_l,n_l2,P_l)=1.0d0/(2.0d0*dble(P_l)-1.0d0)
        end if
        if (isnan(F_new(ind,1,n_l,n_l2,P_l))) then
            print*,'error in generate beliefs'
        end if
        F_new(ind,1:2*P_l-1,n_l,n_l2,P_l)=F_new(ind,1:2*P_l-1,n_l,n_l2,P_l)/sum(F_new(ind,1:2*P_l-1,n_l,n_l2,P_l))
        !print*,sum(F_new(ind,:,n_l,n_l2,P_l))
    end do;end do;end do;end do

!print*,'press key to continue'    
!read*,continue_k
end subroutine
    