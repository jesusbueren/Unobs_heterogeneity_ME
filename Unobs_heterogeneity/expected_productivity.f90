subroutine expected_productivity(params,a,Ef_v,v_l,u_l,re_l)
    use dimensions;use primitives
    implicit none
    double precision,dimension(3),intent(in)::params
    double precision,intent(in)::a
    double precision,dimension(2*P_max-1,3,P_max),intent(out)::Ef_v
    integer,intent(in)::v_l,u_l,re_l
    double precision,dimension(2*P_max,M,3)::Ef
    integer::m_l,k_l,k_l2,n_l,P1_l,P2_l,ind,N,P,p_l
    double precision::theta_p, gamma_p, beta_p,rho_p,cost_p,a_chosen
    
    !Define parameter
    theta_p=params(1)
    beta_p=params(2)
    gamma_p=params(3)
    cost_p=0.0d0!params(4)
    
    !y=theta*( q^beta * a^gamma + re)
    Ef_v=sqrt(-1.0d0)
    !Compute expected productivity and generate the vector of it with the from 2*P-1 form
    do P=1,P_max  
        Ef=0.0d0
        do m_l=1,M
            do p_l=1,2*P-1;do k_l=1,K
                if (theta_p*gamma_p*(q(k_l,1)**beta_p*a**(gamma_p-1.0d0))>cost_p) then
                    a_chosen=a
                else
                    a_chosen=(cost_p/theta_p/gamma_p*q(k_l,1)**beta_p)**(1.0d0/(1.0d0-gamma_p))
                end if
                Ef(p_l,m_l,2)=Ef(p_l,m_l,2)+(theta_p*(q(k_l,1)**beta_p*a_chosen**gamma_p+re_effect(re_l))-cost_p*a_chosen)*PI_k(p_l,k_l,m_l,u_l)
            end do;end do
            do p_l=2,2*P;do k_l=1,K;do k_l2=1,K
                !Ef(2:2*P,m_l,3)=Ef(2:2*P,m_l,3)+theta_p*(q(k_l,1)**(beta_p/(1.0d0-gamma_p))+q(k_l2,1)**(beta_p/(1.0d0-gamma_p)))**(1.0d0-gamma_p)*a**gamma_p &
                !    *PI_k(2:2*P,k_l,m_l,u_l)*PI_k(2:2*P,k_l2,m_l,u_l)
                if (theta_p*gamma_p*((q(k_l,1)+q(k_l2,1))**beta_p*a**(gamma_p-1.0d0))>cost_p) then
                    a_chosen=a
                else
                    a_chosen=(cost_p/theta_p/gamma_p*(q(k_l,1)+q(k_l2,1))**beta_p)**(1.0d0/(1.0d0-gamma_p))
                end if
                Ef(p_l,m_l,3)=Ef(p_l,m_l,3)+(theta_p*(q(k_l,1)+q(k_l2,1))**(beta_p)*a_chosen**gamma_p-cost_p*a_chosen+re_effect(re_l)) &
                    *PI_k(p_l,k_l,m_l,u_l)*PI_k(p_l,k_l2,m_l,u_l)
            end do; end do;end do
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
        !print*,'u_l',u_l
        !print*,'Ef',Ef_v(1:P*2-1,2,P)
    end do
    
    
end subroutine
    