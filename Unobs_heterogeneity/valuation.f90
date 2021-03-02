subroutine valuation(CCP,C,Ef_v,P,V,v_l)
    use dimensions; use primitives
    integer,intent(in)::P
    double precision,dimension((2*P-1)*3,(2*P-1)*3),intent(in)::C
    double precision,dimension(2*P-1,2),intent(in)::CCP
    double precision,dimension(2*P-1,3),intent(in)::Ef_v
    integer,intent(in)::v_l
    double precision,dimension((2*P-1)*3,1),intent(out)::V
    double precision,dimension((2*P-1)*3,1)::U
    double precision,dimension(2*P-1,3)::U_small
    integer::n_l,m_l,ind
    

    U_small(:,3)=0.0d0
    do ind=1,2*P-1; do n_l=1,3
        if (n_l<3) then
            if (CCP(ind,n_l)==1.0d0) then
                U_small(ind,n_l)=(rho*gamma-c_s*PI_s_v(ind,n_l,P,v_l)-c_d*(1.0d0-PI_s_v(ind,n_l,P,v_l)))
            elseif (CCP(ind,n_l)==0.0d0) then
                U_small(ind,n_l)=rho*gamma
            else            
                U_small(ind,n_l)=CCP(ind,n_l)*(rho*gamma-rho*log(CCP(ind,n_l))-c_s*PI_s_v(ind,n_l,P,v_l)-c_d*(1.0d0-PI_s_v(ind,n_l,P,v_l)))+&
                                 (1.0d0-CCP(ind,n_l))*(rho*gamma-rho*log(1.0d0-CCP(ind,n_l)))
            end if
        else
             U_small(ind,n_l)=rho*gamma
        end if
    end do; end do
    U_small=U_small+(1.0d0-tau)*Ef_v
    U_small(:,1)=U_small(:,1)+T_g
       
    do n_l=1,3
        ind=(2*P-1)*(n_l-1)+1
        U(ind:ind+2*P-1-1,1)=U_small(:,n_l)
    end do
    V=matmul(C,U)
   
end subroutine