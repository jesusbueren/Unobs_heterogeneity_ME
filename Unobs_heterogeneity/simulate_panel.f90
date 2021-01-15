subroutine simulate_panel(CCP,n_initial)
    use cadastral_maps; use primitives; use simulation
    implicit none
    double precision,dimension(2*P_max-1,2,P_max,types_a,villages,unobs_types),intent(in)::CCP
    integer,dimension(plots_in_map,villages),intent(inout)::n_initial
    integer::i_l,j_l,t_l,ind,N_all,n_l,P,A,P_l,n_l2,it,v_l,m_l,max_Big_n,s_l,t
    double precision::u_d,u_s,u_f,u_m
    integer,dimension(plots_i)::unobs_types_i
    integer,dimension(1)::seed=456
    double precision,dimension(P_max)::dist
    character::continue_k
    double precision,dimension(types_a)::av_drilling_a,counter_a
    double precision,dimension(2)::av_drilling_t,counter_t
    
    !Call seed number
    call random_seed(PUT=seed)
    
    !Decision to drill
    drilling_it=-9
    counter_a=0.0d0
    counter_t=0.0d0
    
    do s_l=1,simulations
        !Sample unobserved heterogeneity type
        do i_l=1,plots_i
            call RANDOM_NUMBER(u_m)
            if (u_m<UHE_type(1,i_l)) then !these numbers come from hanan's estimation_v3.pdf
                unobs_types_i(i_l)=1
            elseif (u_m<sum(UHE_type(1:2,i_l))) then
                unobs_types_i(i_l)=2
            else
                unobs_types_i(i_l)=3
            end if
        end do
            
        !Store the state for each plot and simulate decision to drill
        do t_l=1,T_sim
            do i_l=1,plots_i
                N_all=-9 !Indicates the number of wells in the adjacency (1 = no wells)
                !Sample the total number of wells 
                ind=1
                do while (N_all==-9)
                    call RANDOM_NUMBER(u_m)
                    if ( (u_m<sum(Pr_N_data(1:ind,t_l,i_l))) .or. ind==max_NFW+1)  then
                        N_all=ind
                    end if
                    ind=ind+1
                end do
                n_l=n_data(t_l,i_l)                 
                P=P_type(i_l)
                A=A_type(i_l)
                v_l=V_type(i_l)
                if (n_l==1) then
                    ind=N_all !position in the state space wrt to the CCP, PI_s_v and ,PI_f_v
                elseif (n_l==2) then
                    ind=N_all-1
                elseif (n_l==3) then
                    ind=N_all-2
                else
                    print*,'error in sim panel'
                end if      
                !Well drilling decision and failures/successes
                if (n_l==1) then !no well
                    call RANDOM_NUMBER(u_d)
                    if (u_d<CCP(ind,n_l,P,A,v_l,unobs_types_i(i_l))) then !decides to drill
                        drilling_it(t_l,i_l,s_l)=1
                    else !decides not to drill
                        drilling_it(t_l,i_l,s_l)=0
                    end if
                elseif (n_l==2) then !one well
                    call RANDOM_NUMBER(u_d)
                    if (u_d<CCP(ind,n_l,P,A,v_l,unobs_types_i(i_l))) then !decides to drill
                        drilling_it(t_l,i_l,s_l)=1
                    else !decides not to drill
                        drilling_it(t_l,i_l,s_l)=0
                    end if 
                elseif(n_l==3) then !two wells
                    drilling_it(t_l,i_l,s_l)=-9
                else
                    print*,'error in simulate panel'
                end if
                
                if (n_l<3) then
                    counter_a(A)=counter_a(A)+1.0d0
                    av_drilling_a(A)=(counter_a(A)-1.0d0)/counter_a(A)*av_drilling_a(A)+1.0d0/counter_a(A)*drilling_it(t_l,i_l,s_l)
                    if (UHE_type(1,i_l)<0.29d0) then
                        t=1
                    else
                        t=2
                    end if
                    counter_t(t)=counter_t(t)+1.0d0
                    av_drilling_t(t)=(counter_t(t)-1.0d0)/counter_t(t)*av_drilling_t(t)+1.0d0/counter_t(t)*drilling_it(t_l,i_l,s_l)
                end if
                
            end do
        end do
    end do
    
    print*,'average drilling across a ',av_drilling_a
    print*,'average drilling across t ',av_drilling_t
    
    OPEN(UNIT=12, FILE="data_sim.txt")
        write(12,*),drilling_it
    close(12)
end subroutine
    