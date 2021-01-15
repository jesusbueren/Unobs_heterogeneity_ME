subroutine expected_productivity(params,a,Ef_v,v_l,u_l)
    use dimensions;use primitives
    implicit none
    double precision,dimension(par),intent(in)::params
    double precision,intent(in)::a
    double precision,dimension(2*P_max-1,3,P_max),intent(out)::Ef_v
    integer,intent(in)::v_l,u_l
    double precision,dimension(2*P_max,M,3)::Ef
    integer::m_l,k_l,k_l2,n_l,P1_l,P2_l,ind,N,P
    double precision::theta_p, gamma_p, beta_p,rho_p
    
    !Define parameter
    theta_p=params(1)
    beta_p=params(2)
    gamma_p=params(3)
    rho_p=1.0d0
    
    !y=theta (beta q^gamma + (1-beta) a^gamma)^(1/(1-gamma))
    Ef_v=sqrt(-1.0d0)
    !Compute expected productivity and generate the vector of it with the from 2*P-1 form
    do P=2,P_max  
        Ef=0.0d0
        do m_l=1,M
            Ef(1:2*P-1,m_l,2)=matmul(PI_k(1:2*P-1,:,m_l,u_l),theta_p*(beta_p*q(:,1)**gamma_p + (1.0d0-beta_p)*a**gamma_p)**(rho_p/gamma_p))
            do k_l=1,K;do k_l2=1,K
                Ef(2:2*P,m_l,3)=Ef(2:2*P,m_l,3)+theta_p*(beta_p*(q(k_l,1)+q(k_l2,1))**gamma_p+ (1.0d0-beta_p)*a**gamma_p)**(rho_p/gamma_p) &
                    *PI_k(2:2*P,k_l,m_l,u_l)*PI_k(2:2*P,k_l2,m_l,u_l)
            end do; end do
        end do 
    
        !Store expected productivity in vector form
        do n_l=1,3
            if (n_l==1) then
                Ef_v(1:P*2-1,n_l,P)=0.0d0
            elseif (n_l==2) then
                Ef_v(1:P*2-1,n_l,P)=PI_m(1,v_l)*Ef(1:2*P-1,1,n_l)+PI_m(2,v_l)*Ef(1:2*P-1,2,n_l) !Compute the expected return across monzoons
            elseif (n_l==3) then
                Ef_v(1:P*2-1,n_l,P)=PI_m(1,v_l)*Ef(2:2*P,1,n_l)+PI_m(2,v_l)*Ef(2:2*P,2,n_l) !Compute the expected return across monzoons
            end if           
        end do
    end do
!Ef(:,:,3)
end subroutine
    