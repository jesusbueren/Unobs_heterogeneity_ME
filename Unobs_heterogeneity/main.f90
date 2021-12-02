program main
use dimensions; use cadastral_maps; use simulation
implicit none
integer,dimension(1)::seed=321
double precision,dimension(par)::params_true,params_MLE
double precision,dimension(2*P_max-1,2*P_max-1,3,3,P_max,villages)::F_true
double precision,dimension(2*P_max-1,2,P_max,types_a,villages,unobs_types)::CCP_true
double precision,dimension(2*P_max-1,3,P_max,types_a,unobs_types)::V_fct
integer,dimension(plots_in_map,villages)::n_dist
double precision,dimension(villages)::mean_N,mean_NPV,mean_budget
double precision::log_likeli
integer::v_l,p_l
character::end_key


!Call seed number
call random_seed(PUT=seed)

!Input primitives of the model
call input_primitives()

!Create a map with plots, neighbors and types
call load_cadastral_maps()

!Load panel data of drilling decisions
call load_estimation_data()

!Uncomment if you want to simulate and check recovery of the simulated parameters.
    !Choose true parameter vector and compute associated equilibrium beliefs and a final distribution of wells in all plots
    !params_true=(/15.66d0,0.7d0,3.2d0/)
    !do v_l=1,villages
    !    print*,'village,',v_l
    !    if (v_l==1) then
    !        CCP_true(:,:,:,:,v_l,:)=0.07d0
    !        n_dist(:,v_l)=1
    !        V_fct=0.0d0
    !    else
    !        CCP_true(:,:,:,:,v_l,:)=CCP_true(:,:,:,:,1,:)
    !        n_dist(:,v_l)=n_dist(:,1)
    !    end if        
    !    call compute_eq_F_CCP(params_true,F_true(:,:,:,:,:,v_l),CCP_true(:,:,:,:,v_l,:),V_fct,n_dist(:,v_l),v_l,mean_N(v_l),mean_NPV(v_l),mean_budget(v_l),Pr_u_X(:,:,:,:,v_l,:))
    !end do
    
    !Using CCPs simulate panel data with drilling decision.
    !OPEN(UNIT=12, FILE="CCP_true.txt")
    !    read(12,*),CCP_true,Pr_u_X
    !close(12)
    !call simulate_panel(CCP_true,n_dist)
    !print*,'end simulate_panel'
    !print*,'Loading panel'
    !OPEN(UNIT=12, FILE="data_sim.txt")
    !    read(12,*),drilling_it
    !close(12)
    
call compute_moments(dble(drilling_it(:,:,1)),"data",moment_own_nxa_data)

print*,'Start estimation'
!Generate a random CCP for computing initial beliefs
CCP_est=sqrt(-1.0d0)
do P_l=1,P_max
    CCP_est(1:2*P_l-1,1:2,P_l,:,:,1)=0.02d0
    CCP_est(1:2*P_l-1,1:2,P_l,:,:,2)=0.07d0
    CCP_est(1:2*P_l-1,1:2,P_l,:,:,3)=0.12d0
end do

call estimation2(params_MLE,log_likeli)


open(unit=12, file=path_results//"parameters.txt")
    read(12,'(<par>f20.12)'),params_mle
close(12)

!call bootstrap_se()
!open(unit=12, file=path_results//"bootstrapped_parameters_85_nem.txt")
!    read(12,*),params_mle
!close(12)
print*,'estimated parameters',params_MLE
!call generate_panel_sample(params_MLE)
call counterfactual_2(params_MLE)
!call counterfactual_1(params_MLE)

!call transitional_dynamics(params_MLE)
!call quantifiying_strategic_interactions(params_MLE)

read*,end_key

end program main