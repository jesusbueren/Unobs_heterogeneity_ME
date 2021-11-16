subroutine update_CCP(V,F,Ef_v,P,CCP,v_l,u_l,a_l)
    use dimensions; use primitives
    implicit none
    integer,intent(in)::P
    double precision,dimension((3*P-2)*3,1),intent(in)::V
    double precision,dimension(3*P-2,3*P-2,4,4),intent(in)::F
    double precision,dimension(3*P-2,4),intent(in)::Ef_v
    integer,intent(in)::v_l,u_l,a_l
    double precision,dimension(3*P-2,3),intent(out)::CCP
    double precision,dimension(3*P-2,4)::V_t
    integer::m_l
    double precision,dimension(3*P-2)::v_00,v_0I,v_10,v_1I,v_20,v_2I
    
    !Vector to tersor
    V_t=reshape(V,(/3*P-2,4/))

        !No well (n=1)
        !!!!!!!!!!!!!!
            !No attempt
            v_00=beta*matmul(F(:,:,1,1),V_t(:,1))
            !Attempt
            v_0I=PI_s_v(1:3*P-2,1,P,v_l)*(-c_s+beta*matmul(F(:,:,1,2),V_t(:,2))) &
                        +(1.0d0-PI_s_v(1:3*P-2,1,P,v_l))*(-c_d+beta*matmul(F(:,:,1,1),V_t(:,1))) 
            
        !!One well (n=2)
        !!!!!!!!!!!!!!!!!
            !No attempt
            v_10=beta*(1.0d0-PI_f_v(1:3*P-2,2,P,v_l,u_l))*matmul(F(:,:,2,2),V_t(:,2)) & !no failure
                +beta*PI_f_v(1:3*P-2,2,P,v_l,u_l)*matmul(F(:,:,2,1),V_t(:,1)) !failure
                
            !Attempt
            v_1I= -PI_s_v(1:3*P-2,2,P,v_l)*c_s-(1.0d0-PI_s_v(1:3*P-2,2,P,v_l))*c_d &
                    +beta*PI_s_v(1:3*P-2,2,P,v_l)*(1.0d0-PI_f_v(1:3*P-2,2,P,v_l,u_l))*matmul(F(:,:,2,3),V_t(:,3)) & !successful & no failure
                    +beta*(PI_s_v(1:3*P-2,2,P,v_l)*PI_f_v(1:3*P-2,2,P,v_l,u_l)+(1.0d0-PI_s_v(1:3*P-2,2,P,v_l))*(1.0d0-PI_f_v(1:3*P-2,2,P,v_l,u_l)))*matmul(F(:,:,2,2),V_t(:,2)) & !successful &  failure or unsuccessful & no failure
                    +beta*(1.0d0-PI_s_v(1:3*P-2,2,P,v_l))*PI_f_v(1:3*P-2,2,P,v_l,u_l)*matmul(F(:,:,2,1),V_t(:,1)) !unsuccessful &  failure
            
        !!Two wells (n=3)
        !!!!!!!!!!!!!!!!!
            !No attempt
            v_20=beta*(1.0d0-PI_f_v(1:3*P-2,3,P,v_l,u_l))**2.0d0*matmul(F(:,:,3,3),V_t(:,3)) & !no failure
                +beta*2.0d0*PI_f_v(1:3*P-2,3,P,v_l,u_l)*(1.0d0-PI_f_v(1:3*P-2,3,P,v_l,u_l))*matmul(F(:,:,3,2),V_t(:,2)) & ! failure of one but not the other
                +beta*PI_f_v(1:3*P-2,3,P,v_l,u_l)**2.0d0*matmul(F(:,:,3,1),V_t(:,1)) ! failure of both
                
            !Attempt
            v_2I= -PI_s_v(1:3*P-2,3,P,v_l)*c_s-(1.0d0-PI_s_v(1:3*P-2,3,P,v_l))*c_d+ &
                beta*(1.0d0-PI_f_v(1:3*P-2,3,P,v_l,u_l))**2.0d0*PI_s_v(1:3*P-2,3,P,v_l)*matmul(F(:,:,3,4),V_t(:,4)) & !no failure & success
                +beta*(2.0d0*PI_f_v(1:3*P-2,3,P,v_l,u_l)*(1.0d0-PI_f_v(1:3*P-2,3,P,v_l,u_l))*PI_s_v(1:3*P-2,3,P,v_l)+(1.0d0-PI_f_v(1:3*P-2,3,P,v_l,u_l))**2.0d0*(1.0d0-PI_s_v(1:3*P-2,3,P,v_l)))*matmul(F(:,:,3,3),V_t(:,3)) & ! failure of one & success or failure of non & unsuccess
                +beta*(PI_f_v(1:3*P-2,3,P,v_l,u_l)**2.0d0*PI_s_v(1:3*P-2,3,P,v_l)+2.0d0*PI_f_v(1:3*P-2,3,P,v_l,u_l)*(1.0d0-PI_f_v(1:3*P-2,3,P,v_l,u_l))*(1.0d0-PI_s_v(1:3*P-2,3,P,v_l)))*matmul(F(:,:,3,2),V_t(:,2)) & ! failure of both & succesful or failure of one and unsucessful
                +beta*(PI_f_v(1:3*P-2,3,P,v_l,u_l)**2.0d0*(1.0d0-PI_s_v(1:3*P-2,3,P,v_l)))*matmul(F(:,:,3,1),V_t(:,1)) ! failure of two & unsucessful
                    

    
    CCP(:,1)=1.0d0/(1.0d0+exp(v_00/rho(1)-v_0I/rho(1)))
    CCP(:,2)=1.0d0/(1.0d0+exp(v_10/rho(2)-v_1I/rho(2)))
    CCP(:,3)=1.0d0/(1.0d0+exp(v_20/rho(2)-v_2I/rho(2)))
    
end subroutine