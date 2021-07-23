!Estimate parameters
subroutine estimation(params_MLE,log_likeli)
    use dimensions; use primitives; use simulation
    implicit none
    double precision,dimension(par),intent(out)::params_MLE
    double precision,intent(out)::log_likeli
    double precision::log_L,ftol
    double precision,dimension(par+1,par)::p_g
    double precision,dimension(par+1)::y
    integer::iter,p_l,a_l,P_l2,v_l,ind,n_l,u_l
    interface 
        function log_likelihood(params_MLE)
            use dimensions
            double precision,dimension(par),intent(in)::params_MLE
            double precision::log_likelihood
        end function log_likelihood
    end interface
    integer,dimension(plots_in_map,villages)::n_initial_all
    double precision,dimension(2*P_max-1,2,P_max,types_a,villages,unobs_types)::CCP_old,CCP_mid
    double precision,dimension(2*P_max-1,3,P_max,types_a,villages,unobs_types)::V_fct
    double precision,dimension(2*P_max-1,3,P_max,types_a,villages,unobs_types)::Ef_v !Ef_v: expected productivity
    double precision::dist
    integer::it
    integer(8),dimension(2*P_max-1,3,3,P_max,villages)::iterations_all
    double precision,dimension(par,par)::xi
    integer,dimension(1)::seed_c
    double precision, dimension(villages)::mean_N,mean_NPV,mean_budget
    character::pause_k
    
    it=1
    iterations_all=0.0d0
    
    CCP_mid=CCP_est
    !Generate beliefs consitent with CCP
    F_est=1.0d0
    Ef_v=0.0d0
   print*,'Generating beliefs'
    !Generate an initial well endowment: everyone has zero wells
1   n_initial_all(1:plots_in_map,:)=1    
    call random_seed(GET=seed_c)
    !$OMP PARALLEL default(shared) 
    !$OMP  DO
    do v_l=1,villages
        print*,'village ',v_l,' out of ',villages
        call generate_beliefs(CCP_mid(:,:,:,:,v_l,:),V_fct(:,:,:,:,v_l,:),Ef_v(:,:,:,:,v_l,:),n_initial_all(:,v_l),F_est(:,:,:,:,:,v_l),v_l,iterations_all(:,:,:,:,v_l),mean_N(v_l),mean_NPV(v_l),mean_budget(v_l),Pr_u_X(:,:,:,:,v_l,:))
    end do
    !$OMP END DO  
    !$OMP END PARALLEL      
    
    call random_seed(PUT=seed_c)
    !Fixing beliefs, estimate parameter
    !print*,'Initial Conditions'
    
    if (it==1) then
        p_g(1,:)=(/6.11d0,0.5d0,0.5d0,5.0d0/) !(/3.03d0,0.10d0,0.51d0,15.04d0/)
    end if
    
    do p_l=2,par+1
        p_g(p_l,:)=p_g(1,:)
        p_g(p_l,p_l-1)=p_g(1,p_l-1)*0.8d0
    end do
        
    !Change parameters to the (-Inf;Inf) real line
    do p_l=1,par+1
        p_g(p_l,1)=log(p_g(p_l,1))
        p_g(p_l,2:3)=log(p_g(p_l,2:3)/(1.0d0-p_g(p_l,2:3)))
        p_g(p_l,4)=log(p_g(p_l,4))
        y(p_l)=log_likelihood(p_g(p_l,:))
        !print*,'press key to continue'
        !print*,'press key to continue'
        !read*,pause_k 
    end do 
    !print*,'likelihood_ini',y(1)
    
    !Optimization of parameters given beliefs
    print*,'game against nature'
    ftol=1.0d-6
    call amoeba(p_g,y,ftol,log_likelihood,iter)
    print*,'likelihood amoeba',y(1)
    p_g(:,1)=exp(p_g(:,1))
    p_g(:,2:3)=1.0d0/(1.0d0 + exp(-p_g(:,2:3))) 
    p_g(:,4)=exp(p_g(:,4))
    print*,' parameters amoeba',p_g(1,:)
    !Change parameters to the (-Inf;Inf) real line
    p_g(1,1)=log(p_g(1,1))
    p_g(1,2:3)=log(p_g(1,2:3)/(1.0d0-p_g(1,2:3)))
    p_g(1,4)=log(p_g(1,4))
    xi=0.0d0
    do p_l=1,par
        xi(p_l,p_l)=1.0d0
    end do
    ftol=1.0d-4
    !call powell(p_g(1,:),xi,ftol,iter,y(1))
    log_likeli=y(1)
    p_g(:,1)=exp(p_g(:,1))
    p_g(:,2:3)=1.0d0/(1.0d0 + exp(-p_g(:,2:3))) 
    p_g(:,4)=exp(p_g(:,4))
    !print*,'likelihood powell',y(1)
    !print*,'parameter powell',p_g(1,:)
    
    
    !Compute CCP to check convergence
    params_MLE=p_g(1,:)
    CCP_old=CCP_est
    rho=params_MLE(4)
    
    do v_l=1,villages
        do u_l=1,unobs_types;do a_l=1,types_a
            call expected_productivity(params_MLE(1:3),area(a_l),Ef_v(:,:,:,a_l,v_l,u_l),v_l,u_l)
        end do;end do
        do P_l=1,P_max; do a_l=1,types_a; do u_l=1,unobs_types 
            call policy_fct_it(Ef_v(1:2*P_l-1,:,P_l,a_l,v_l,u_l)&
                                ,F_est(1:2*P_l-1,1:2*P_l-1,:,:,P_l,v_l) &
                                ,P_l &
                                ,CCP_old(1:2*P_l-1,:,P_l,a_l,v_l,u_l),CCP_est(1:2*P_l-1,:,P_l,a_l,v_l,u_l),v_l,u_l &
                                 ,V_fct(1:2*P_l-1,:,P_l,a_l,v_l,u_l))
            !call value_fct_it(Ef_v(1:2*P_l-1,:,P_l,a_l,v_l,u_l)&
            !                ,F_est(1:2*P_l-1,1:2*P_l-1,:,:,P_l,v_l) &
            !                ,P_l &
            !                ,CCP_est(1:2*P_l-1,:,P_l,a_l,v_l,u_l),v_l,u_l &
            !                ,V_fct(1:2*P_l-1,:,P_l,a_l,v_l,u_l)) 
        end do; end do;end do
    end do
    
    !print*,'est CCP',CCP_est(1:2*P_max-1,:,P_max,:,:,:)

    
    dist=0.0
    do P_l=1,P_max; do n_l=1,2;do ind=1,2*P_l-1; do v_l=1,villages
        dist=dist+dble(sum(iterations_all(ind,n_l,1:3,P_l,v_l)))/dble(sum(iterations_all(:,1:2,1:3,:,:)))*sum(abs(CCP_old(ind,n_l,P_l,:,v_l,:)-CCP_est(ind,n_l,P_l,:,v_l,:)))/dble(types_a)/dble(unobs_types)
    end do;end do; end do;end do
    print*,'dist',dist
    
    !New guess of the ccp is half way through
    CCP_mid=CCP_est*0.5d0+CCP_old*0.5d0
    if (dist>5.0d-4) then
        it=it+1
        go to 1
    end if
    
    if (dist<0.0) then
        print*,'error in precision of iterations'   
    end if

   
end subroutine
    
function log_likelihood(params_MLE)
    use dimensions; use simulation; use cadastral_maps; use primitives
    implicit none
    double precision,dimension(par),intent(in)::params_MLE
    double precision,dimension(par)::params
    double precision,dimension(2*P_max-1,2,P_max,types_a,villages,unobs_types)::CCP
    double precision,dimension(2*P_max-1,3,P_max,types_a,villages,unobs_types)::V_fct
    integer::i_l,t_l,type_l,a_l,p_l,v_l,ind,u_l,j_l,s_l,t,missing_x1
    double precision::log_likelihood
    double precision,dimension(2*P_max-1,3,P_max,types_a,unobs_types)::Ef_v !Ef_v: expected productivity
    double precision,dimension(unobs_types)::likelihood_i,likelihood_it
    character::end_k
    double precision,dimension(unobs_types)::av_CCP_uhe
    double precision,dimension(T_sim,plots_i)::av_CCP_it
    double precision,dimension(plots_i)::likelihood_aux
    character::pause_k
    
    
    params(1)=exp(params_MLE(1))
    params(2:3)=1.0d0/(1.0d0 + exp(-params_MLE(2:3))) 
    params(4)=exp(params_MLE(4))
    
    rho=exp(params_MLE(4))
    print*,' parameters',params
    
    log_likelihood=0.0d0
    missing_x1=0
    
    do a_l=1,types_a; do v_l=1,villages;do u_l=1,unobs_types
        call expected_productivity(params(1:3),area(a_l),Ef_v(:,:,:,a_l,u_l),v_l,u_l)
    end do; end do;end do
    

    !$OMP PARALLEL default(shared) 
    !$OMP  DO
    do P_l=1,P_max; do a_l=1,types_a ; do u_l=1,unobs_types;do v_l=1,villages
        call policy_fct_it(Ef_v(1:2*P_l-1,:,P_l,a_l,u_l)&
                        ,F_est(1:2*P_l-1,1:2*P_l-1,:,:,P_l,v_l) &
                        ,P_l &
                        ,CCP_est(1:2*P_l-1,:,P_l,a_l,v_l,u_l),CCP(1:2*P_l-1,:,P_l,a_l,v_l,u_l),v_l,u_l,V_fct(1:2*P_l-1,:,P_l,a_l,v_l,u_l))
        !call value_fct_it(Ef_v(1:2*P_l-1,:,P_l,a_l,u_l)&
        !                    ,F_est(1:2*P_l-1,1:2*P_l-1,:,:,P_l,v_l) &
        !                    ,P_l &
        !                    ,CCP(1:2*P_l-1,:,P_l,a_l,v_l,u_l),v_l,u_l &
        !                    ,V_fct(1:2*P_l-1,:,P_l,a_l,v_l,u_l)) 
    end do; end do;end do; end do
    !$OMP END DO  
    !$OMP END PARALLEL

    do s_l=1,simulations
    do i_l=1,plots_i;
        if (impute_i(i_l)==0) then
            UHE_type_model(:,i_l)=0.0d0
            if (P_type(i_l)>1) then !more than one neighbor
                likelihood_i=1.0d0
                do t_l=1,T_sim
                    likelihood_it=0.0d0
                    av_CCP_uhe=0.0d0
                    !if (sum(Pr_N_data(n_data(t_l,i_l):min(max_NFW+1,2*(P_type(i_l)-1)+n_data(t_l,i_l)),t_l,i_l))<0.98d0) then
                    !    print*,'error in estimation'
                    !    print*,Pr_N_data(n_data(t_l,i_l):min(max_NFW+1,2*(P_type(i_l)-1)+n_data(t_l,i_l)),t_l,i_l)
                    !    read*,pause_k
                    !end if
                
                    !Set distribution of unobserved heterogeneity
                    if (t_l==1) then
                        do j_l=n_data(t_l,i_l),min(max_NFW+1,2*(P_type(i_l)-1)+n_data(t_l,i_l))
                            if (n_data(t_l,i_l)==1) then
                                ind=j_l !position in the state space wrt to the CCP, PI_s_v and ,PI_f_v
                            elseif (n_data(t_l,i_l)==2) then
                                ind=j_l-1
                            elseif (n_data(t_l,i_l)==3) then
                                ind=j_l-2
                            else
                                print*,'error in estimation'
                            end if 
                            UHE_type_model(:,i_l)=UHE_type_model(:,i_l)+Pr_u_X(ind,n_data(t_l,i_l),P_type(i_l),A_type(i_l),V_type(i_l),:)*Pr_N_data(j_l,t_l,i_l)
                        end do
                    end if
                
                    do j_l=n_data(t_l,i_l),min(max_NFW+1,2*(P_type(i_l)-1)+n_data(t_l,i_l))
                        if (n_data(t_l,i_l)==1) then
                            ind=j_l !position in the state space wrt to the CCP, PI_s_v and ,PI_f_v
                        elseif (n_data(t_l,i_l)==2) then
                            ind=j_l-1
                        elseif (n_data(t_l,i_l)==3) then
                            ind=j_l-2
                        else
                            print*,'error in estimation'
                        end if
                        if (drilling_it(t_l,i_l,s_l)==1) then
                            likelihood_it=likelihood_it+CCP(ind,n_data(t_l,i_l),P_type(i_l),A_type(i_l),V_type(i_l),:)*Pr_N_data(j_l,t_l,i_l)
                            av_CCP_uhe=av_CCP_uhe+CCP(ind,n_data(t_l,i_l),P_type(i_l),A_type(i_l),V_type(i_l),:)*Pr_N_data(j_l,t_l,i_l)
                        elseif (drilling_it(t_l,i_l,s_l)==0) then
                            likelihood_it=likelihood_it+(1.0d0-CCP(ind,n_data(t_l,i_l),P_type(i_l),A_type(i_l),V_type(i_l),:))*Pr_N_data(j_l,t_l,i_l)
                            av_CCP_uhe=av_CCP_uhe+CCP(ind,n_data(t_l,i_l),P_type(i_l),A_type(i_l),V_type(i_l),:)*Pr_N_data(j_l,t_l,i_l)
                        else
                            likelihood_it=likelihood_it+Pr_N_data(j_l,t_l,i_l) !so that if drilling_it is missing likelihood_it is one
                        end if
                    end do
                    likelihood_i=likelihood_i*likelihood_it
                    if (isnan(sum(likelihood_i))) then
                        print*,'pb in likelihood',i_l,ind,n_data(t_l,i_l),P_type(i_l),A_type(i_l),V_type(i_l),CCP(ind,n_data(t_l,i_l),P_type(i_l),A_type(i_l),V_type(i_l),2),drilling_it(t_l,i_l,s_l)
                        read*,end_k
                    end if
                    if (drilling_it(t_l,i_l,s_l)==1 .or. drilling_it(t_l,i_l,s_l)==0) then
                        if (UHE_type(1,i_l)==-9.0d0) then
                            print*,n_data(t_l,i_l)
                        end if
                        av_CCP_it(t_l,i_l)=sum(av_CCP_uhe*UHE_type(:,i_l))
                    else
                        av_CCP_it(t_l,i_l)=-9.0d0
                    end if
                end do;
                if (UHE_type_model(1,i_l)==-9.0d0) then
                    UHE_type_model(:,i_l)=1.0d0/dble(unobs_types)
                end if
                if (sum(UHE_type_model(:,i_l))/=0.0d0)then
                    log_likelihood=log_likelihood+log(sum(likelihood_i*UHE_type_model(:,i_l)))
                else
                    missing_x1=missing_x1+1
                end if
                likelihood_aux(i_l)=log(sum(likelihood_i*UHE_type_model(:,i_l)))
                !if (log(sum(likelihood_i*UHE_type_model(:,i_l)))==-1.0/0.0) then
                !    print*,CCP(ind,n_data(t_l,i_l),P_type(i_l),A_type(i_l),V_type(i_l),:)
                !    print*,'paused'
                !    read*,pause_k
                !end if
            end if
        end if
    end do
    end do
        
    log_likelihood=-log_likelihood
    
    if (bootstrap==0) then
        call compute_moments(av_CCP_it,"modl")
    end if
    
    print*,'likelihood',log_likelihood
    print*,'missing_x1',missing_x1
    !print*,'paused'
    !read*,pause_k
end function