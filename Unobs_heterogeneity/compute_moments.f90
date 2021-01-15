subroutine compute_moments(data_in,string_name)
    use simulation
    implicit none
    double precision,dimension(T_sim,plots_i),intent(in)::data_in
    CHARACTER (LEN=4),intent(in) :: string_name
    double precision,dimension(max_NFW+1)::counter_N,moment_N
    double precision,dimension(types_a,2)::counter_own_nxa,moment_own_nxa
    double precision,dimension(unobs_types,2)::counter_uhe,moment_uhe
    double precision,dimension(P_max)::counter_P,moment_P
    integer::i_l,t_l
    
    !Initialize to zero
    counter_N=0.0d0
    moment_N=0.0d0
    counter_own_nxa=0.0d0
    moment_own_nxa=0.0d0
    counter_uhe=0.0d0
    moment_uhe=0.0d0
    counter_P=0.0d0
    moment_P=0.0d0
    
    do i_l=1,plots_i;do t_l=1,T_sim
        if (data_in(t_l,i_l)/=-9.0d0) then
            !Moments across number of functioning wells in adjacency
            counter_N(modal_N(t_l,i_l))=counter_N(modal_N(t_l,i_l))+1.0d0
            moment_N(modal_N(t_l,i_l))=(counter_N(modal_N(t_l,i_l))-1.0)/counter_N(modal_N(t_l,i_l))*moment_N(modal_N(t_l,i_l))&
                                        +1.0d0/counter_N(modal_N(t_l,i_l))*data_in(t_l,i_l)
            !Moments across number of owned functioning wells and owned wells
            counter_own_nxa(A_type(i_l),n_data(t_l,i_l))=counter_own_nxa(A_type(i_l),n_data(t_l,i_l))+1.0d0
            moment_own_nxa(A_type(i_l),n_data(t_l,i_l))=(counter_own_nxa(A_type(i_l),n_data(t_l,i_l))-1.0)/counter_own_nxa(A_type(i_l),n_data(t_l,i_l))*moment_own_nxa(A_type(i_l),n_data(t_l,i_l))&
                                            +1.0d0/counter_own_nxa(A_type(i_l),n_data(t_l,i_l))*data_in(t_l,i_l)
            !Moments across number of unobserved heterogeneity types
            counter_uhe(modal_UHE_type(i_l),n_data(t_l,i_l))=counter_uhe(modal_UHE_type(i_l),n_data(t_l,i_l))+1.0d0
            moment_uhe(modal_UHE_type(i_l),n_data(t_l,i_l))=(counter_uhe(modal_UHE_type(i_l),n_data(t_l,i_l))-1.0)/counter_uhe(modal_UHE_type(i_l),n_data(t_l,i_l))*moment_uhe(modal_UHE_type(i_l),n_data(t_l,i_l))&
                                            +1.0d0/counter_uhe(modal_UHE_type(i_l),n_data(t_l,i_l))*data_in(t_l,i_l)
            !Moments across number of plots
            counter_P(P_type(i_l))=counter_P(P_type(i_l))+1.0d0
            moment_P(P_type(i_l))=(counter_P(P_type(i_l))-1.0)/counter_P(P_type(i_l))*moment_P(P_type(i_l))&
                                    +1.0d0/counter_P(P_type(i_l))*data_in(t_l,i_l)  
        end if
        
    end do;end do
    
    
    OPEN(UNIT=12, FILE=path_results//string_name//"_N.txt")
        write(12,*),moment_N
    close(12)
    
    OPEN(UNIT=12, FILE=path_results//string_name//"_own_nxa.txt")
        write(12,*),moment_own_nxa
    close(12)
    
    OPEN(UNIT=12, FILE=path_results//string_name//"_uhe.txt")
        write(12,*),moment_uhe
    close(12)
    
    OPEN(UNIT=12, FILE=path_results//string_name//"_P.txt")
        write(12,*),moment_P
    close(12)
    
    
end subroutine
    
    