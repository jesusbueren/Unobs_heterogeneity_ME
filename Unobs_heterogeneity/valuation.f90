subroutine valuation(CCP,C,Ef_v,P,V,v_l,u_l,a_l,CCP_aux)
    use dimensions; use primitives
    integer,intent(in)::P
    double precision,dimension((3*P-2)*4,(3*P-2)*4),intent(in)::C
    double precision,dimension(3*P-2,3),intent(in)::CCP
    double precision,dimension(3*P-2),intent(in)::CCP_aux
    double precision,dimension(3*P-2,4),intent(in)::Ef_v
    integer,intent(in)::v_l,u_l,a_l
    double precision,dimension((3*P-2)*4,1),intent(out)::V
    double precision,dimension((3*P-2)*4,1)::U
    double precision,dimension(3*P-2,4)::U_small
    integer::n_l,m_l,ind
    

    !U_small=0.0d0
    !do ind=1,3*P-2; do n_l=1,3
    !    if (n_l==1) then
    !        if (CCP(ind,n_l)==1.0d0)then
    !            U_small(ind,n_l)=v_nod+rho(n_l)*gamma-c_s*PI_s_v(ind,n_l,P,v_l)-c_d*(1.0d0-PI_s_v(ind,n_l,P,v_l))              
    !        elseif (CCP(ind,n_l)==0.0d0)then
    !            U_small(ind,n_l)=v_nod+rho(n_l)*gamma
    !        else
    !            U_small(ind,n_l)=v_nod+CCP(ind,n_l)*(rho(n_l)*gamma-rho(n_l)*log(CCP(ind,n_l))-c_s*PI_s_v(ind,n_l,P,v_l)-c_d*(1.0d0-PI_s_v(ind,n_l,P,v_l)))+&
    !                            (1.0d0-CCP(ind,n_l))*(rho(n_l)*gamma-rho(n_l)*log(1.0d0-CCP(ind,n_l)))
    !        end if
    !    elseif (n_l==2) then
    !        if (CCP(ind,n_l)==1.0d0)then
    !            U_small(ind,n_l)=rho(n_l)*gamma-c_s*PI_s_v(ind,n_l,P,v_l)-c_d*(1.0d0-PI_s_v(ind,n_l,P,v_l))              
    !        elseif (CCP(ind,n_l)==0.0d0)then
    !            U_small(ind,n_l)=rho(n_l)*gamma
    !        else
    !            U_small(ind,n_l)=CCP(ind,n_l)*(rho(n_l)*gamma-rho(n_l)*log(CCP(ind,n_l))-c_s*PI_s_v(ind,n_l,P,v_l)-c_d*(1.0d0-PI_s_v(ind,n_l,P,v_l)))+&
    !                            (1.0d0-CCP(ind,n_l))*(rho(n_l)*gamma-rho(n_l)*log(1.0d0-CCP(ind,n_l)))
    !        end if
    !    else
    !        if (CCP_aux(ind)==1.0d0 )then
    !            U_small(ind,n_l)=rho(2)*gamma-c_s*PI_s_v(ind,2,P,v_l)-c_d*(1.0d0-PI_s_v(ind,2,P,v_l))
    !        elseif (CCP_aux(ind)==0.0d0 )then
    !            U_small(ind,n_l)=rho(2)*gamma
    !        else
    !             U_small(ind,n_l)=CCP_aux(ind)*(rho(2)*gamma-rho(2)*log(CCP_aux(ind))-c_s*PI_s_v(ind,2,P,v_l)-c_d*(1.0d0-PI_s_v(ind,2,P,v_l)))+&
    !                                (1.0d0-CCP_aux(ind))*(rho(2)*gamma-rho(2)*log(1.0d0-CCP_aux(ind))) !CCP(ind,2)*(rho(2)*gamma-rho(2)*log(CCP(ind,2))-c_s*PI_s_v(ind,2,P,v_l)-c_d*(1.0d0-PI_s_v(ind,2,P,v_l)))+&
    !                                !(1.0d0-CCP(ind,2))*(rho(2)*gamma-rho(2)*log(1.0d0-CCP(ind,2))) 
    !         end if
    !    end if
    !end do; end do
    
    U_small=0.0d0
    do ind=1,3*P-2; do n_l=1,3
        if (n_l==1) then
            if (CCP(ind,n_l)==1.0d0)then
                U_small(ind,n_l)=-c_s*PI_s_v(ind,n_l,P,v_l)-c_d*(1.0d0-PI_s_v(ind,n_l,P,v_l))              
            elseif (CCP(ind,n_l)==0.0d0)then
                U_small(ind,n_l)=v_nod
            else
                U_small(ind,n_l)=CCP(ind,n_l)*(-c_s*PI_s_v(ind,n_l,P,v_l)-c_d*(1.0d0-PI_s_v(ind,n_l,P,v_l)))
            end if
        elseif (n_l==2 .or. n_l==3) then
            if (CCP(ind,n_l)==1.0d0)then
                U_small(ind,n_l)=-c_s*PI_s_v(ind,n_l,P,v_l)-c_d*(1.0d0-PI_s_v(ind,n_l,P,v_l))              
            elseif (CCP(ind,n_l)==0.0d0)then
                U_small(ind,n_l)=0.0d0
            else
                U_small(ind,n_l)=CCP(ind,n_l)*(-c_s*PI_s_v(ind,n_l,P,v_l)-c_d*(1.0d0-PI_s_v(ind,n_l,P,v_l)))
            end if
        else
            if (CCP_aux(ind)==1.0d0 )then
                U_small(ind,n_l)=0.0d0
            elseif (CCP_aux(ind)==0.0d0 )then
                U_small(ind,n_l)=0.0d0
            else
                 U_small(ind,n_l)=0.0d0
             end if
        end if
    end do; end do
    
    U_small=U_small+Ef_v
    !Tax on having a well
    if (social==1) then
        U_small(:,2)=U_small(:,2)-c_e
        U_small(:,3)=U_small(:,3)-2.0d0*c_e
        U_small(:,4)=U_small(:,4)-3.0d0*c_e
    end if
    
    
    do n_l=1,3
        ind=(3*P-2)*(n_l-1)+1
        U(ind:ind+3*P-2-1,1)=U_small(:,n_l)
    end do
    V=matmul(C,U)
   
    
end subroutine