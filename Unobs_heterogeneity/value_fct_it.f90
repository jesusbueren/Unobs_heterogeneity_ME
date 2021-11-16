subroutine value_fct_it(Ef_v,F,P,CCP,v_l,u_l,V_new)
    use dimensions;use primitives
    implicit none
    integer,intent(in)::P
    double precision,dimension(3*P-2,3),intent(in)::Ef_v
    double precision,dimension(3*P-2,3*P-2,4,4),intent(in)::F
    integer,intent(in)::v_l,u_l
    double precision,dimension(3*P-2,4)::V_old
    double precision,dimension((3*P-2)*4)::Vec_old,Vec_new
    double precision::crit,dist
    double precision,dimension(3*P-2,3),intent(out)::CCP
    double precision,dimension(3*P-2,4),intent(out)::V_new
    character::pause_k

    !Define intial guess of value function
    V_old=0.0d0
    dist=1.0d0
    crit=1.0d-10
        
    do while (dist>crit)
        call one_step_value_fct_it(Ef_v,F,P,CCP,v_l,u_l,V_old,V_new)
                
        !Check contraction
        !!!!!!!!!!!!!!!!!!
            Vec_old=reshape(V_old,(/(3*P-2)*4/))
            Vec_new=reshape(V_new,(/(3*P-2)*4/))
            dist=maxval(abs(Vec_old-Vec_new))
            if (dist==1.0d0/0.0d0)then
                print*,'error in vfl it'
                print*,'Vec_old',Vec_old
                print*,'Vec_new',Vec_new
                read*,pause_k
            end if
            if (dist>crit)then
                V_old=V_new
            end if
    end do   
end subroutine
    
        
subroutine one_step_value_fct_it(Ef_v,F,P,CCP,v_l,u_l,V_old,V_new)
    use dimensions;use primitives
    implicit none
    integer,intent(in)::P
    double precision,dimension(3*P-2,3),intent(in)::Ef_v
    double precision,dimension(3*P-2,3*P-2,4,4),intent(in)::F
    integer,intent(in)::v_l,u_l
    double precision,dimension(3*P-2,4),intent(in)::V_old
    double precision,dimension(3*P-2,4),intent(out)::V_new
    double precision,dimension((3*P-2)*4)::Vec_old,Vec_new
    double precision,dimension(3*P-2)::v_00,v_0I,v_10,v_1I,v_20,v_2I,CCP_aux
    double precision,dimension(3*P-2,3),intent(out)::CCP
    character::pause_k
    
    V_new=0.0d0
    CCP=0.0d0
    CCP_aux=1.0d0/(1.0d0+exp(-(-PI_s_v(1:3*P-3,3,P,v_l)*c_s-(1.0d0-PI_s_v(1:3*P-3,3,P,v_l))*c_d)/rho(2)))
        
            !No well (n=1)
            !!!!!!!!!!!!!!
                !No attempt
                v_00(1:3*P-2)=T_g+beta*matmul(F(1:3*P-2,1:3*P-2,1,1),V_old(1:3*P-2,1))
                !Attempt               
                v_0I(1:3*P-2)=T_g+PI_s_v(1:3*P-2,1,P,v_l)*(-c_s+beta*(matmul(F(1:3*P-2,1:3*P-2,1,2),V_old(1:3*P-2,2)))) & !success
                           +(1.0d0-PI_s_v(1:3*P-2,1,P,v_l))*(-c_d+beta*(matmul(F(1:3*P-2,1:3*P-2,1,1),V_old(1:3*P-2,1)))) !Failure
                !Value function
                V_new(1:3*P-2,1)=rho(1)*log(1.0d0+exp((v_0I(1:3*P-2)-v_00(1:3*P-2))/rho(1)))+v_00(1:3*P-2)+rho(1)*gamma

            !!One well (n=2)
            !!!!!!!!!!!!!!!!!
                !No attempt
                v_10(1:3*P-2)=T_g+Ef_v(1:3*P-2,2)-tau &
                    +beta*((1.0d0-PI_f_v(1:3*P-2,2,P,v_l,u_l))*matmul(F(1:3*P-2,1:3*P-2,2,2),V_old(1:3*P-2,2))) & !No failure
                    +beta*(PI_f_v(1:3*P-2,2,P,v_l,u_l)*matmul(F(1:3*P-2,1:3*P-2,2,1),V_old(1:3*P-2,1))) !Failure
    
                !Attempt
                v_1I(1:3*P-2)=T_g+Ef_v(1:3*P-2,2)-tau &
                        -PI_s_v(1:3*P-2,2,P,v_l)*c_s-(1.0d0-PI_s_v(1:3*P-2,2,P,v_l))*c_d &
                        +beta*(PI_s_v(1:3*P-2,2,P,v_l)*(1.0d0-PI_f_v(1:3*P-2,2,P,v_l,u_l))*matmul(F(1:3*P-2,1:3*P-2,2,3),V_old(1:3*P-2,3))) & !Success in the new and no failure of the old
                        +beta*(PI_s_v(1:3*P-2,2,P,v_l)*PI_f_v(1:3*P-2,2,P,v_l,u_l)+(1.0d0-PI_s_v(1:3*P-2,2,P,v_l))*(1.0d0-PI_f_v(1:3*P-2,2,P,v_l,u_l)))*matmul(F(1:3*P-2,1:3*P-2,2,2),V_old(1:3*P-2,2)) & ! Success and failure of the old or failure of the new but no faile of the old
                        +beta*((1.0d0-PI_s_v(1:3*P-2,2,P,v_l))*PI_f_v(1:3*P-2,2,P,v_l,u_l)*matmul(F(1:3*P-2,1:3*P-2,2,1),V_old(1:3*P-2,1))) !failure in both
                
                !Value function
                V_new(1:3*P-2,2)=rho(2)*log(1.0d0+exp((v_1I(1:3*P-2)-v_10(1:3*P-2))/rho(2)))+v_10(1:3*P-2)+rho(2)*gamma
            
            !!Two wells (n=3)
            !!!!!!!!!!!!!!!!!!
                !No attempt
                v_20(1:3*P-2)=T_g+Ef_v(1:3*P-2,3)-2.0d0*tau &
                                + beta*(1.0d0-PI_f_v(1:3*P-2,3,P,v_l,u_l))**2.0d0*matmul(F(1:3*P-2,1:3*P-2,3,3),V_old(1:3*P-2,3)) & !none fails
                                + 2.0d0*beta*(1.0d0-PI_f_v(1:3*P-2,3,P,v_l,u_l))*PI_f_v(1:3*P-2,3,P,v_l,u_l)*matmul(F(1:3*P-2,1:3*P-2,3,2),V_old(1:3*P-2,2)) & !one fails
                                + beta*PI_f_v(1:3*P-2,3,P,v_l,u_l)**2.0d0*matmul(F(1:3*P-2,1:3*P-2,3,1),V_old(1:3*P-2,1))  !both fail
    
                !Attempt
                v_20(1:3*P-2)=T_g+Ef_v(1:3*P-2,3)-2.0d0*tau &
                                -PI_s_v(1:3*P-2,3,P,v_l)*c_s-(1.0d0-PI_s_v(1:3*P-2,3,P,v_l))*c_d &
                                + beta*(1.0d0-PI_f_v(1:3*P-2,3,P,v_l,u_l))**2.0d0*PI_s_v(1:3*P-2,3,P,v_l)*matmul(F(1:3*P-2,1:3*P-2,3,4),V_old(1:3*P-2,4)) & !none fails & success
                                + 2.0d0*beta*(1.0d0-PI_f_v(1:3*P-2,3,P,v_l,u_l))*PI_f_v(1:3*P-2,3,P,v_l,u_l)*PI_s_v(1:3*P-2,3,P,v_l)*matmul(F(1:3*P-2,1:3*P-2,3,3),V_old(1:3*P-2,3)) & !one fails & success
                                + beta*PI_f_v(1:3*P-2,3,P,v_l,u_l)**2.0d0*PI_f_v(1:3*P-2,3,P,v_l,u_l)*PI_s_v(1:3*P-2,3,P,v_l)*matmul(F(1:3*P-2,1:3*P-2,3,2),V_old(1:3*P-2,2)) &  !both fail & success
                                + beta*(1.0d0-PI_f_v(1:3*P-2,3,P,v_l,u_l))**2.0d0*(1.0d0-PI_s_v(1:3*P-2,3,P,v_l))*matmul(F(1:3*P-2,1:3*P-2,3,3),V_old(1:3*P-2,3)) & !none fails & no success
                                + 2.0d0*beta*(1.0d0-PI_f_v(1:3*P-2,3,P,v_l,u_l))*PI_f_v(1:3*P-2,3,P,v_l,u_l)*(1.0d0-PI_s_v(1:3*P-2,3,P,v_l))*matmul(F(1:3*P-2,1:3*P-2,3,2),V_old(1:3*P-2,2)) & !one fails & no success
                                + beta*PI_f_v(1:3*P-2,3,P,v_l,u_l)**2.0d0*PI_f_v(1:3*P-2,3,P,v_l,u_l)*(1.0d0-PI_s_v(1:3*P-2,3,P,v_l))*matmul(F(1:3*P-2,1:3*P-2,3,1),V_old(1:3*P-2,1))  !both fail & no success
                
                !Value function
                V_new(1:3*P-2,3)=rho(2)*log(1.0d0+exp((v_2I(1:3*P-2)-v_20(1:3*P-2))/rho(2)))+v_20(1:3*P-2)+rho(2)*gamma
                
                
            !!Three wells (n=4)
            !!!!!!!!!!!!!!!!!!
                CCP(1:3*P-2,1)=1.0d0/(1.0d0+exp(v_00(1:3*P-2)/rho(2)-v_0I(1:3*P-2)/rho(2)))
                CCP(1:3*P-2,2)=1.0d0/(1.0d0+exp(v_10(1:3*P-2)/rho(2)-v_1I(1:3*P-2)/rho(2)))
                CCP(1:3*P-2,3)=1.0d0/(1.0d0+exp(v_10(1:3*P-2)/rho(2)-v_1I(1:3*P-2)/rho(2)))
                
                V_new(1:3*P-2,4)=T_g+Ef_v(1:3*P-2,3)-3.0d0*tau &
                                +CCP_aux(1:3*P-2)*(rho(2)*gamma-rho(2)*log(CCP_aux(1:3*P-2))-c_s*PI_s_v(1:3*P-2,3,P,v_l)-c_d*(1.0d0-PI_s_v(1:3*P-2,3,P,v_l)))&
                                +(1.0d0-CCP_aux(1:3*P-2))*(rho(2)*gamma-rho(2)*log(1.0d0-CCP_aux(1:3*P-2))) &
                                + beta*(1.0d0-PI_f_v(1:3*P-2,3,P,v_l,u_l))**3.0d0*matmul(F(1:3*P-2,1:3*P-2,4,4),V_old(1:3*P-2,4)) & !none fails
                                + beta*3.0d0*PI_f_v(1:3*P-2,3,P,v_l,u_l)*(1.0d0-PI_f_v(1:3*P-2,3,P,v_l,u_l))**2.0d0*matmul(F(1:3*P-2,1:3*P-2,4,3),V_old(1:3*P-2,3)) & !one fails
                                + beta*3.0d0*PI_f_v(1:3*P-2,3,P,v_l,u_l)**2.0d0*(1.0d0-PI_f_v(1:3*P-2,3,P,v_l,u_l))*matmul(F(1:3*P-2,1:3*P-2,4,2),V_old(1:3*P-2,2)) & !two fail
                                + beta*PI_f_v(1:3*P-2,3,P,v_l,u_l)**3.0d0*matmul(F(1:3*P-2,1:3*P-2,4,1),V_old(1:3*P-2,1))   !three fail
    
    
    if (isnan(sum(CCP))) then
        print*,'error in vfi',v_00(1:3*P-2)
        read*,pause_k
    end if
    
    

end subroutine