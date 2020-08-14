% =========================================================================
% Part of the Matlab code to accompany the paper
% 'Solving heterogeneous agent models in discrete time with
% many idiosyncratic states by perturbation methods', CEPR Discussion Paper
% DP13071
% by Christian Bayer and Ralph Luetticke
% http://www.erc-hump-uni-bonn.net/
% =========================================================================

%% Initialize workspace and load directories

addpath(genpath('functions'))
addpath(genpath('latex'))

%% Switch options
casename='SS_BASELINE_HANC_JPEG';
printIRFs     = false;
mpar.overrideEigen = false;

%% Produce matrices to reduce state-space
% disp('Reduce State Space and Compute System for Steady State')
mainskript_statereduc

F = @(a,b,c,d)Fsys(a,b,c,d,Xss,Yss,Gamma_state,indexMUdct,DC1,DC2,Copula,par,mpar,grid,meshes,P_H,aggrshock,oc,os);

[Fss,LHS,RHS,Distr] = F(State,State_m,Contr,Contr_m);


%% Solve RE via Schmitt-Grohe-Uribe Form
% disp('Take Numerical Derivatives and Solve RE via Schmitt-Grohe-Uribe Form')

[hx,gx,F1,F2,F3,F4,par] = SGU_solver(F, mpar.numstates,mpar.numcontrols,oc,mpar,par,p);

%%
eta0=zeros(mpar.numstates,1);
eta0(end)=1;

[gxx_gss,hxx_hss,gss_L,hss_L,F11,F12,F13,F14,F21,F22,F23,F24,F31,F32,F33,F34,F41,F42,F43,F44,par] = ...
       SGU_solver_2nd_L(F,F1,F2,F3,F4,gx,hx,mpar.numstates,mpar.numcontrols,oc,mpar,par,grid,eta0);
   
   %%
    for i=1:mpar.numstates
    gxx(:,:,i) = gxx_gss(:,(i-1)*mpar.numstates+i:i*mpar.numstates+i-1);
    hxx(:,:,i) = hxx_hss(:,(i-1)*mpar.numstates+i:i*mpar.numstates+i-1);
  end
   
  gss = gxx_gss(:,end);
  hss = hxx_hss(:,end);
  
