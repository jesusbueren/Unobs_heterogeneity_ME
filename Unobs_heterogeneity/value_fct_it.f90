subroutine value_fct_it(Ef_v,F,P,CCP,v_l,u_l,V_new)
    use dimensions;use primitives
    implicit none
    integer,intent(in)::P
    double precision,dimension(2*P-1,3),intent(in)::Ef_v
    double precision,dimension(2*P-1,2*P-1,3,3),intent(in)::F
    integer,intent(in)::v_l,u_l
    double precision,dimension(2*P-1,3)::V_old
    double precision,dimension((2*P-1)*3)::Vec_old,Vec_new
    double precision,dimension(2*P-1)::v_00,v_0I,v_10,v_1I
    double precision::crit,dist
    integer::n_l,m_l
    double precision,dimension(2*P-1,2),intent(out)::CCP
    double precision,dimension(2*P-1,3),intent(out)::V_new
    character::pause_k

    
    
    !Define intial guess of value function
    V_old=0.0d0
    dist=1.0d0
    crit=1.0d-15
    
    
    do while (dist>crit)
            !No well (n=1)
            !!!!!!!!!!!!!!
                !No attempt
                v_00(1:2*P-1)=T_g+beta*matmul(F(1:2*P-1,1:2*P-1,1,1),V_old(1:2*P-1,1))
                !Attempt               
                v_0I(1:2*P-1)=T_g+PI_s_v(1:2*P-1,1,P,v_l)*(-c_s+beta*(matmul(F(1:2*P-1,1:2*P-1,1,2),V_old(1:2*P-1,2)))) & !success
                           +(1.0d0-PI_s_v(1:2*P-1,1,P,v_l))*(-c_d+beta*(matmul(F(1:2*P-1,1:2*P-1,1,1),V_old(1:2*P-1,1)))) !Failure
                !Value function
                V_new(1:2*P-1,1)=rho*log(1.0d0+exp((v_0I(1:2*P-1)-v_00(1:2*P-1))/rho))+v_00(1:2*P-1)

            !!One well (n=2)
            !!!!!!!!!!!!!!!!!
                !No attempt
                v_10(1:2*P-1)=T_g+Ef_v(1:2*P-1,2) &
                    +beta*((1.0d0-PI_f_v(1:2*P-1,2,P,v_l,u_l))*matmul(F(1:2*P-1,1:2*P-1,2,2),V_old(1:2*P-1,2))) & !No failure
                    +beta*(PI_f_v(1:2*P-1,2,P,v_l,u_l)*matmul(F(1:2*P-1,1:2*P-1,2,1),V_old(1:2*P-1,1))) !Failure

                !Attempt
                v_1I(1:2*P-1)=T_g+Ef_v(1:2*P-1,2) &
                        -PI_s_v(1:2*P-1,2,P,v_l)*c_s-(1.0d0-PI_s_v(1:2*P-1,2,P,v_l))*c_d &
                        +beta*(PI_s_v(1:2*P-1,2,P,v_l)*(1.0d0-PI_f_v(1:2*P-1,2,P,v_l,u_l))*matmul(F(1:2*P-1,1:2*P-1,2,3),V_old(1:2*P-1,3))) & !Success in the new and no failure of the old
                        +beta*(PI_s_v(1:2*P-1,2,P,v_l)*PI_f_v(1:2*P-1,2,P,v_l,u_l)+(1.0d0-PI_s_v(1:2*P-1,2,P,v_l))*(1.0d0-PI_f_v(1:2*P-1,2,P,v_l,u_l)))*matmul(F(1:2*P-1,1:2*P-1,2,2),V_old(1:2*P-1,2)) & ! Success and failure of the old or failure of the new but no faile of the old
                        +beta*((1.0d0-PI_s_v(1:2*P-1,2,P,v_l))*PI_f_v(1:2*P-1,2,P,v_l,u_l)*matmul(F(1:2*P-1,1:2*P-1,2,1),V_old(1:2*P-1,1))) !failure in both
                
                !Value function
                V_new(1:2*P-1,2)=rho*log(1.0d0+exp((v_1I(1:2*P-1)-v_10(1:2*P-1))/rho))+v_10(1:2*P-1)
            
            !!Two wells (n=3)
            !!!!!!!!!!!!!!!!!!
                V_new(1:2*P-1,3)=T_g+Ef_v(1:2*P-1,3) &
                                + beta*(1.0d0-PI_f_v(1:2*P-1,3,P,v_l,u_l))**2.0d0*matmul(F(1:2*P-1,1:2*P-1,3,3),V_old(1:2*P-1,3)) & !none fails
                                + 2.0d0*beta*(1.0d0-PI_f_v(1:2*P-1,3,P,v_l,u_l))*PI_f_v(1:2*P-1,3,P,v_l,u_l)*matmul(F(1:2*P-1,1:2*P-1,3,2),V_old(1:2*P-1,2)) & !one fails
                                + beta*PI_f_v(1:2*P-1,3,P,v_l,u_l)**2.0d0*matmul(F(1:2*P-1,1:2*P-1,3,1),V_old(1:2*P-1,1))  !both fail


                
        !Check contraction
        !!!!!!!!!!!!!!!!!!
            Vec_old=reshape(V_old,(/(2*P-1)*3/))
            Vec_new=reshape(V_new,(/(2*P-1)*3/))
            dist=maxval(abs(Vec_old-Vec_new))
            if (dist==1.0d0/0.0d0)then
                print*,Vec_old
                read*,pause_k
            end if
            if (dist>crit)then
                V_old=V_new
            end if
    end do

    CCP(1:2*P-1,1)=1.0d0/(1.0d0+exp(v_00(1:2*P-1)/rho-v_0I(1:2*P-1)/rho))
    CCP(1:2*P-1,2)=1.0d0/(1.0d0+exp(v_10(1:2*P-1)/rho-v_1I(1:2*P-1)/rho))
    
    if (isnan(sum(CCP))) then
        print*,'error in vfi',v_00(1:2*P-1)
        read*,pause_k
    end if
    
    
    V_new=V_new+(rho*gamma)/(1.0d0-beta)
    
    

end subroutine