subroutine generate_C(CCP,F,P,C,v_l,u_l)
    use dimensions; use primitives
    implicit none
    integer,intent(in)::P
    double precision,dimension(3*P-2,3),intent(in)::CCP
    double precision,dimension(3*P-2,3*P-2,4,4),intent(in)::F
    integer,intent(in)::v_l,u_l
    integer::m_l,n_l,m_l2,n_l2,ind,ind2,i,j
    double precision,dimension(1,3*P-2)::ones
    double precision,dimension(3*P-2,3*P-2,4,4)::F_big
    double precision,dimension((3*P-2)*4,(3*P-2)*4)::super_F
    double precision,dimension((3*P-2)*4,(3*P-2)*4),intent(out)::C
    double precision,dimension((3*P-2)*4,(3*P-2)*4)::C_old,Iden
    double precision,dimension((3*P-2)*4*(3*P-2)*4)::C_old_v,C_v
    double precision::dist,crit
    
    ones=1.0d0
    F_big(:,:,1,1)=matmul(CCP(:,1:1)*(1.0d0-PI_s_v(1:3*P-2,1:1,P,v_l))+(1.0d0-CCP(:,1:1)),ones) &
                            *F(:,:,1,1)
    F_big(:,:,1,2)=matmul(CCP(:,1:1)*PI_s_v(1:3*P-2,1:1,P,v_l),ones) &
                            *F(:,:,1,2)
    F_big(:,:,1,3)=0.0d0
    F_big(:,:,2,1)=matmul((CCP(:,2:2)*(1.0d0-PI_s_v(1:3*P-2,2:2,P,v_l))+(1.0d0-CCP(:,2:2)))*PI_f_v(1:3*P-2,2:2,P,v_l,u_l),ones) &
                            *F(:,:,2,1)
    F_big(:,:,2,2)=matmul((CCP(:,2:2)*(1.0d0-PI_s_v(1:3*P-2,2:2,P,v_l))+(1.0d0-CCP(:,2:2)))*(1.0d0-PI_f_v(1:3*P-2,2:2,P,v_l,u_l))+&
                            CCP(:,2:2)*PI_s_v(1:3*P-2,2:2,P,v_l)*PI_f_v(1:3*P-2,2:2,P,v_l,u_l),ones) &
                            *F(:,:,2,2)
    F_big(:,:,2,3)=matmul(CCP(:,2:2)*PI_s_v(1:3*P-2,2:2,P,v_l)*(1.0d0-PI_f_v(1:3*P-2,2:2,P,v_l,u_l)),ones) &
                            *F(:,:,2,3)
    F_big(:,:,3,1)=matmul((1.0d0-CCP(:,3:3)+CCP(:,2:2)*(1.0d0-PI_s_v(1:3*P-2,3:3,P,v_l)))*PI_f_v(1:3*P-2,3:3,P,v_l,u_l)**2.0d0,ones) &
                            *F(:,:,3,1)
    F_big(:,:,3,2)=2.0d0*matmul((1.0d0-CCP(:,3:3)+CCP(:,3:3)*(1.0d0-PI_s_v(1:3*P-2,3:3,P,v_l)))*PI_f_v(1:3*P-2,3:3,P,v_l,u_l)*(1.0d0-PI_f_v(1:3*P-2,3:3,P,v_l,u_l))+ &
                                CCP(:,3:3)*PI_s_v(1:3*P-2,3:3,P,v_l)*PI_f_v(1:3*P-2,3:3,P,v_l,u_l)**2.0d0,ones) &
                            *F(:,:,3,2)
    F_big(:,:,3,3)=matmul((1.0d0-CCP(:,3:3))*(1.0d0-PI_f_v(1:3*P-2,3:3,P,v_l,u_l))*(1.0d0-PI_f_v(1:3*P-2,3:3,P,v_l,u_l)) + &
                            2.0d0*CCP(:,3:3)*PI_s_v(1:3*P-2,3:3,P,v_l)*PI_f_v(1:3*P-2,3:3,P,v_l,u_l)*(1.0d0-PI_f_v(1:3*P-2,3:3,P,v_l,u_l)),ones) &
                            *F(:,:,3,3)
    F_big(:,:,4,1)=matmul(PI_f_v(1:3*P-2,3:3,P,v_l,u_l)**3.0d0,ones) &
                            *F(:,:,4,1)
    F_big(:,:,4,2)=matmul(3.0d0*PI_f_v(1:3*P-2,3:3,P,v_l,u_l)**2.0d0*(1.0d0-PI_f_v(1:3*P-2,3:3,P,v_l,u_l)),ones) &
                            *F(:,:,4,2)
    F_big(:,:,4,3)=matmul(3.0d0*PI_f_v(1:3*P-2,3:3,P,v_l,u_l)*(1.0d0-PI_f_v(1:3*P-2,3:3,P,v_l,u_l))**2.0d0,ones) &
                            *F(:,:,4,3)
    F_big(:,:,4,4)=matmul((1.0d0-PI_f_v(1:3*P-2,3:3,P,v_l,u_l))**3.0d0,ones) &
                            *F(:,:,4,4)

    do n_l=1,3
        ind=(3*P-2)*(n_l-1)+1
        do n_l2=1,3;
            ind2=(3*P-2)*(n_l2-1)+1
            super_F(ind:ind+(3*P-2)-1,ind2:ind2+(3*P-2)-1)=F_big(:,:,n_l,n_l2)
        end do
    end do

    Iden=0.0d0
    do i=1,(3*P-2)*4; do j=1,(3*P-2)*4
        if (i==j)then
            Iden(i,j)=1.0d0
        end if
    end do; end do
    
    C_old=Iden
    dist=1.0d0
    crit=1.0d-9
    do while (dist>crit)
        C=Iden+beta*matmul(super_F,C_old)
        C_old_v=reshape(C_old,(/(3*P-2)*4*(3*P-2)*4/))
        C_v=reshape(C,(/(3*P-2)*3*(3*P-2)*4/))
        dist= maxval(sqrt((C_old_v-C_v)**2.0d0))
        if (dist>crit)then
            C_old=C
        end if        
    end do
    
end subroutine