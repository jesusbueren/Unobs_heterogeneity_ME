subroutine valuation(CCP,C,Ef_v,P,V,v_l,u_l,a_l)
    use dimensions; use primitives
    integer,intent(in)::P
    double precision,dimension((2*P-1)*3,(2*P-1)*3),intent(in)::C
    double precision,dimension(2*P-1,2),intent(in)::CCP
    double precision,dimension(2*P-1,3),intent(in)::Ef_v
    integer,intent(in)::v_l,u_l,a_l
    double precision,dimension((2*P-1)*3,1),intent(out)::V
    double precision,dimension((2*P-1)*3,1)::U
    double precision,dimension(2*P-1,3)::U_small
    integer::n_l,m_l,ind
    

    U_small=0.0d0
    do ind=1,2*P-1; do n_l=1,3
        if (n_l==1) then
            if (CCP(ind,n_l)==1.0d0)then
                U_small(ind,n_l)=v_nod+rho(n_l)*gamma-c_s*PI_s_v(ind,n_l,P,v_l)-c_d*(1.0d0-PI_s_v(ind,n_l,P,v_l))              
            elseif (CCP(ind,n_l)==0.0d0)then
                U_small(ind,n_l)=v_nod+rho(n_l)*gamma
            else
                U_small(ind,n_l)=v_nod+CCP(ind,n_l)*(rho(n_l)*gamma-rho(n_l)*log(CCP(ind,n_l))-c_s*PI_s_v(ind,n_l,P,v_l)-c_d*(1.0d0-PI_s_v(ind,n_l,P,v_l)))+&
                                (1.0d0-CCP(ind,n_l))*(rho(n_l)*gamma-rho(n_l)*log(1.0d0-CCP(ind,n_l)))
            end if
        elseif (n_l==2) then
            if (CCP(ind,n_l)==1.0d0)then
                U_small(ind,n_l)=rho(n_l)*gamma-c_s*PI_s_v(ind,n_l,P,v_l)-c_d*(1.0d0-PI_s_v(ind,n_l,P,v_l))              
            elseif (CCP(ind,n_l)==0.0d0)then
                U_small(ind,n_l)=rho(n_l)*gamma
            else
                U_small(ind,n_l)=CCP(ind,n_l)*(rho(n_l)*gamma-rho(n_l)*log(CCP(ind,n_l))-c_s*PI_s_v(ind,n_l,P,v_l)-c_d*(1.0d0-PI_s_v(ind,n_l,P,v_l)))+&
                                (1.0d0-CCP(ind,n_l))*(rho(n_l)*gamma-rho(n_l)*log(1.0d0-CCP(ind,n_l)))
            end if
        else
            !if (CCP(ind,2)==1.0d0 )then
            !    U_small(ind,n_l)=rho(2)*gamma-c_s*PI_s_v(ind,2,P,v_l)-c_d*(1.0d0-PI_s_v(ind,2,P,v_l))
            !elseif (CCP(ind,2)==0.0d0 )then
            !    U_small(ind,n_l)=rho(2)*gamma
            !else
                 U_small(ind,n_l)=smthg(ind,P,a_l,v_l,u_l)*(rho(2)*gamma-rho(2)*log(smthg(ind,P,a_l,v_l,u_l))-c_s*PI_s_v(ind,2,P,v_l)-c_d*(1.0d0-PI_s_v(ind,2,P,v_l)))+&
                                    (1.0d0-smthg(ind,P,a_l,v_l,u_l))*(rho(2)*gamma-rho(2)*log(1.0d0-smthg(ind,P,a_l,v_l,u_l))) !CCP(ind,2)*(rho(2)*gamma-rho(2)*log(CCP(ind,2))-c_s*PI_s_v(ind,2,P,v_l)-c_d*(1.0d0-PI_s_v(ind,2,P,v_l)))+&
                                    !(1.0d0-CCP(ind,2))*(rho(2)*gamma-rho(2)*log(1.0d0-CCP(ind,2))) 
             !end if
        end if
    end do; end do
    U_small=U_small+Ef_v+T_g
    !Tax on having a well
    U_small(:,2)=U_small(:,2)-tau
    U_small(:,3)=U_small(:,3)-2.0d0*tau
    
    
    
    do n_l=1,3
        ind=(2*P-1)*(n_l-1)+1
        U(ind:ind+2*P-1-1,1)=U_small(:,n_l)
    end do
    V=matmul(C,U)
   
    
end subroutine