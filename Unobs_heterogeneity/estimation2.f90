!Estimate parameters
subroutine estimation2(params_MLE,log_likeli)
    use dimensions; use primitives; use simulation
    implicit none
    double precision,dimension(par),intent(out)::params_MLE
    double precision,intent(out)::log_likeli
    double precision::log_L,ftol
    double precision,dimension(par+1,par)::p_g
    double precision,dimension(par+1)::y
    integer::iter,p_l,a_l,P_l2,v_l,ind,n_l,u_l
    interface 
        function log_likelihood2(params_MLE)
            use dimensions
            double precision,dimension(par),intent(in)::params_MLE
            double precision::log_likelihood2
        end function log_likelihood2
    end interface
    double precision,dimension(2*P_max-1,2,P_max,types_a,villages,unobs_types)::CCP_old,CCP_mid
    double precision,dimension(2*P_max-1,3,P_max,types_a,villages,unobs_types)::V_fct,V_social
    double precision,dimension(2*P_max-1,3,P_max,types_a,villages,unobs_types)::Ef_v !Ef_v: expected productivity
    double precision::dist
    integer::it
    integer(8),dimension(2*P_max-1,3,3,P_max,villages,unobs_types)::iterations_all
    double precision,dimension(par,par)::xi
    integer,dimension(1)::seed_c
    double precision, dimension(villages)::mean_N,mean_NPV,mean_budget
    double precision,dimension(villages)::village_fe
    character::pause_k
    
      
    call random_seed(PUT=seed_c)
    !Fixing beliefs, estimate parameter
    !print*,'Initial Conditions'
    

    p_g(1,1:4)=(/12.33d0,0.24d0,0.35d0,10.5d0/) !(/16.36d0,0.33d0,0.14d0,8.0d0/)
    do p_l=2,par+1
        p_g(p_l,:)=p_g(1,:)
        p_g(p_l,p_l-1)=p_g(1,p_l-1)*0.8d0
    end do

        
    !Change parameters to the (-Inf;Inf) real line
    do p_l=1,par+1
        p_g(p_l,1)=log(p_g(p_l,1))
        p_g(p_l,2:3)=log(p_g(p_l,2:3)/(1.0d0-p_g(p_l,2:3)))
        p_g(p_l,4)=log(p_g(p_l,4))
        y(p_l)=log_likelihood2(p_g(p_l,:))  
        read*,pause_k
    end do 

    !print*,'likelihood_ini',y(1)
        
    ftol=1.0d-7
    max_mle=99999999.0d0
    call amoeba(p_g,y,ftol,log_likelihood2,iter)
    
    log_likeli=y(1)
    p_g(:,1)=exp(p_g(:,1))
    p_g(:,2:3)=1.0d0/(1.0d0+exp(-p_g(:,2:3)))
    p_g(:,4)=exp(p_g(:,4))

    
    params_MLE=p_g(1,:)
    print*,'likelihood amoeba',y(1)
    
    
    if (dist<0.0) then
        print*,'error in precision of iterations'   
    end if

   
end subroutine
    
function log_likelihood2(params_MLE)
    use dimensions; use simulation; use cadastral_maps; use primitives
    implicit none
    double precision,dimension(par),intent(in)::params_MLE
    double precision,dimension(par)::params
    double precision,dimension(2*P_max-1,2,P_max,types_a,villages,unobs_types)::CCP
    double precision,dimension(2*P_max-1,3,P_max,types_a,villages,unobs_types)::V_fct,V_social
    integer::i_l,t_l,type_l,a_l,p_l,v_l,ind,u_l,j_l,s_l,t,missing_x1,j_l2,ind2
    double precision::log_likelihood2,pr_non_zombie_II
    double precision,dimension(2*P_max-1,3,P_max,types_a,villages,unobs_types)::Ef_v !Ef_v: expected productivity Ef_v(1,:,1,1,1,1)
    double precision,dimension(unobs_types)::likelihood_i,likelihood_it,P_N2_N1,P_BigN2_BigN1
    character::end_k
    double precision,dimension(T_sim,plots_i,unobs_types)::av_CCP_uhe
    double precision,dimension(T_sim,plots_i)::av_CCP_it
    double precision,dimension(plots_i)::likelihood_aux
    character::pause_k
    double precision,dimension(types_a,2)::moment_own_nxa_model
    double precision,dimension(villages)::village_fe
    double precision,dimension(COV,1)::X
    double precision,dimension(2*P_max-1,2*P_max-1,3,3,P_max,villages,unobs_types)::F
    integer,dimension(plots_in_map,villages)::n_dist
    double precision,dimension(villages)::mean_N,social_output,private_output
    double precision,dimension(3,3,2*P_max-1,2,P_max,types_a,villages,unobs_types)::joint_pr
    double precision,dimension(unobs_types)::expost_types
    
    params(1)=exp(params_MLE(1))
    params(2:3)=1.0d0/(1.0d0 + exp(-params_MLE(2:3))) 
    params(4)=exp(params_MLE(4))
    rho=params(4)

    !params(5)=1.0d0/(1.0d0 + exp(-params_MLE(5))) 
    !params(6:18)=exp(params_MLE(6:18))
    
    pr_non_zombie_II=1.0d0
    village_fe=1.0d0
    !village_fe(2:villages)=params(6:18)
    


    print*,' parameters',params
    
    log_likelihood2=0.0d0
    missing_x1=0
    

    CCP=0.07d0
    n_dist=1
    !$OMP PARALLEL default(shared)
    !$OMP  DO
    do v_l=1,1!villages
        print*,'village',v_l        
        call compute_eq_F_CCP(params,F(:,:,:,:,:,v_l,:),CCP(:,:,:,:,v_l,:),V_fct(:,:,:,:,v_l,:),V_social(:,:,:,:,v_l,:),n_dist(:,v_l),v_l,mean_N(v_l),social_output(v_l),private_output(v_l),joint_pr(:,:,:,:,:,:,v_l,:))
    end do
    !$OMP END DO  
    !$OMP END PARALLEL 
     open(unit=12, file=path_results//"init_cond.txt")   
        write(12,*),joint_pr,CCP
    close(12)
    !do v_l=1,villages
    !    CCP(:,:,:,:,v_l,:)=CCP(:,:,:,:,1,:)
    !    joint_pr(:,:,:,:,:,:,v_l,:)=joint_pr(:,:,:,:,:,:,1,:)
    !end do


    

    
    do s_l=1,simulations
    do i_l=1,plots_i;
        X(:,1)=(/1.0d0,N_bar(i_l),dble(A_type(i_l)-1),dble(P_type(i_l))/)
        if (impute_i(i_l)==0) then
            UHE_type_model(:,i_l)=0.0d0
            if (P_type(i_l)>1) then !more than one neighbor
                likelihood_i=1.0d0
                do t_l=1,T_sim
                    likelihood_it=0.0d0
                    av_CCP_uhe(t_l,i_l,:)=0.0d0

                    do j_l=n_data(t_l,i_l),min(max_NFW+1,2*(P_type(i_l)-1)+n_data(t_l,i_l))
                        !position in the state space wrt to the CCP, PI_s_v and ,PI_f_v
                        if (n_data(t_l,i_l)==1) then
                            ind=j_l 
                        elseif (n_data(t_l,i_l)==2) then
                            ind=j_l-1
                        elseif (n_data(t_l,i_l)==3) then
                            ind=j_l-2
                        else
                            print*,'error in estimation'
                        end if
                        if (drilling_it(t_l,i_l,s_l)==1 .and. n_data(t_l,i_l)<3) then 
                            likelihood_it=likelihood_it+CCP(ind,n_data(t_l,i_l),P_type(i_l),A_type(i_l),V_type(i_l),:)* &
                                            sum(reshape(joint_pr(n_data(t_l,i_l),:,ind,:,P_type(i_l),A_type(i_l),V_type(i_l),:),(/3*2,unobs_types/)),1)
                            av_CCP_uhe(t_l,i_l,:)=av_CCP_uhe(t_l,i_l,:)+CCP(ind,n_data(t_l,i_l),P_type(i_l),A_type(i_l),V_type(i_l),:)* &
                                            sum(reshape(joint_pr(n_data(t_l,i_l),:,ind,:,P_type(i_l),A_type(i_l),V_type(i_l),:),(/3*2,unobs_types/)),1)
                        elseif (drilling_it(t_l,i_l,s_l)==0 .and. n_data(t_l,i_l)<3) then
                            likelihood_it=likelihood_it+(1.0d0-CCP(ind,n_data(t_l,i_l),P_type(i_l),A_type(i_l),V_type(i_l),:))*&
                                            sum(reshape(joint_pr(n_data(t_l,i_l),:,ind,:,P_type(i_l),A_type(i_l),V_type(i_l),:),(/3*2,unobs_types/)),1)
                            av_CCP_uhe(t_l,i_l,:)=av_CCP_uhe(t_l,i_l,:)+CCP(ind,n_data(t_l,i_l),P_type(i_l),A_type(i_l),V_type(i_l),:)* &
                                            sum(reshape(joint_pr(n_data(t_l,i_l),:,ind,:,P_type(i_l),A_type(i_l),V_type(i_l),:),(/3*2,unobs_types/)),1)
                        !elseif (drilling_it(t_l,i_l,s_l)==0 .and. n_data(t_l,i_l)==3) then
                        !    likelihood_it=likelihood_it+(1.0d0-CCP_aux(ind,V_type(i_l)))*Pr_N_data(j_l,t_l,i_l)
                        else
                            likelihood_it=1.0d0 !so that if drilling_it is missing likelihood_it is one
                        end if
                    end do
                    
                    if (t_l>1) then
                        if (n_data(t_l-1,i_l)<3) then
                            likelihood_it=likelihood_it*sum(joint_pr(n_data(t_l-1,i_l),n_data(t_l,i_l),:,drilling_it(t_l-1,i_l,s_l)+1,P_type(i_l),A_type(i_l),V_type(i_l),:),1)
                        else
                            likelihood_it=likelihood_it*sum(joint_pr(n_data(t_l-1,i_l),n_data(t_l,i_l),:,1,P_type(i_l),A_type(i_l),V_type(i_l),:),1)
                        end if
                    else
                        likelihood_it=likelihood_it*sum(reshape(joint_pr(:,n_data(t_l,i_l),:,:,P_type(i_l),A_type(i_l),V_type(i_l),:),(/3*(2*P_max-1)*2,3/)),1)
                    end if   

                    likelihood_i=likelihood_i*likelihood_it

                    if (likelihood_i(1)==0.0d0)then
                        print*,''
                    end if


                    !if (isnan(sum(likelihood_i))) then
                    !    print*,'pb in likelihood',i_l,ind,n_data(t_l,i_l),P_type(i_l),A_type(i_l),V_type(i_l),CCP(ind,n_data(t_l,i_l),P_type(i_l),A_type(i_l),V_type(i_l),2),drilling_it(t_l,i_l,s_l)
                    !    read*,end_k
                    !end if
                    !if (sum(likelihood_i)==0.0d0) then
                    !    print*,'pb in likelihood',i_l,ind,n_data(t_l,i_l),P_type(i_l),A_type(i_l),V_type(i_l),CCP(ind,n_data(t_l,i_l),P_type(i_l),A_type(i_l),V_type(i_l),2),drilling_it(t_l,i_l,s_l)
                    !    read*,end_k
                    !end if
                        
                end do;
                

                !Model 5
                if (sum(likelihood_i)/=0.0d0 )then 
                    if (can_be_zombie_i(i_l)==0) then
                        log_likelihood2=log_likelihood2+log(sum(likelihood_i*pr_unobs_t)*pr_non_zombie_II)  
                    else
                        log_likelihood2=log_likelihood2+log(sum(likelihood_i*pr_unobs_t)*pr_non_zombie_II+(1.0d0-pr_non_zombie_II)) 
                    end if
                else
                    missing_x1=missing_x1+1
                end if


                do t_l=1,T_sim
                    if ((drilling_it(t_l,i_l,s_l)==1 .or. drilling_it(t_l,i_l,s_l)==0) .and. sum(likelihood_i)/=0.0d0) then
                        if (can_be_zombie_i(i_l)==0) then
                            av_CCP_it(t_l,i_l)=sum(av_CCP_uhe(t_l,i_l,:)*(likelihood_i/sum(likelihood_i))) 
                        else
                            av_CCP_it(t_l,i_l)=sum(av_CCP_uhe(t_l,i_l,:)*(likelihood_i*pr_unobs_t*pr_non_zombie_II/sum(likelihood_i*pr_unobs_t*pr_non_zombie_II+(1-pr_non_zombie_II))))
                        end if    
                    else
                        av_CCP_it(t_l,i_l)=-9.0d0
                    end if
                end do
                        
            end if
        end if
    end do
    end do
        !close(12)

    log_likelihood2=-log_likelihood2

    
    !GMM
    !log_likelihood=sum(((moment_own_nxa_model-moment_own_nxa_data))**2)
    if (bootstrap==0 .and. log_likelihood2<max_mle) then
        open(unit=12, file=path_results//"parameters.txt",status='replace')
            write(12,'(<par>f20.12,f20.12)'),params,log_likelihood2
        close(12)
        call compute_moments(av_CCP_it,"modl",moment_own_nxa_model)
        max_mle=log_likelihood2
    end if
        
    
    print*,'likelihood',log_likelihood2
    print*,'missing_x1',missing_x1
    !print*,'paused'
    !read*,pause_k
end function