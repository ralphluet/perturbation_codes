% =========================================================================
% Part of the Matlab code to accompany the paper
% 'Solving heterogeneous agent models in discrete time with
% many idiosyncratic states by perturbation methods', CEPR Discussion Paper
% DP13071
% by Christian Bayer and Ralph Luetticke
% http://www.erc-hump-uni-bonn.net/
% =========================================================================
%
% Permission is hereby granted, free of charge, to any person obtaining a
% copy of this software and associated documentation files (the "Software"),
% to deal in the Software without restriction, including without limitation
% the rights to use, copy, modify, merge, publish, distribute, sublicense,
% and/or sell copies of the Software, and to permit persons to whom the
% Software is furnished to do so, subject to the following conditions:
%
% The most recent version or successor of the above-mentioned paper is
% properly cited in all work and publications that benefit from insights
% derived from the Software.
% The above copyright notice and this permission notice shall be included
% in all copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
% THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
% FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
% DEALINGS IN THE SOFTWARE.
%
% =========================================================================

%% Initialize workspace and load directories
clear
clc
close all
restoredefaultpath

Computername='METHOD' 
starttime=clock;
casenamemaster='SS_BASELINE_HANC_PERT';
addpath(genpath('shared_functions'))
oldpath = path;


p = gcp('nocreate');
pool='local'; % select machine

if isempty(p)
    c = parcluster;
    parpool(pool,c.NumWorkers)
    p = gcp('nocreate');
end

%% Select aggregate shock
aggrshock           = 'TFP';
par.rhoS            = 0.75;     % Persistence 
par.sigmaS          = 0.007;    % STD 
mpar.ns             = 3;        % number of TFP states
mpar.T(1,1)         = 1500;     % number of periods
mpar.T(2,1)         = 500;      % burn-in phase
mpar.maxlag         = 40;       % number of periods IRF
mpar.T_KS           = mpar.T(1,1);

mpar.nK             = 3;        % number of grid points for KS


%% Draw random seed for simulation across models
pr_s_master=rand(RandStream('mcg16807', 'Seed',20180621),1,mpar.T(1,1));

[P_S, grid.s] = rouwen(par.rhoS, 0, par.sigmaS/sqrt(1-par.rhoS^2), mpar.ns);
grid.s=exp(grid.s');
PS_S=cumsum(P_S,2);

smaster(1)        =ceil(mpar.ns/2);
for t=1:mpar.T(1,1) -1  
    smaster(t+1) = min(mpar.ns,sum(PS_S(smaster(t),:)<pr_s_master(t))+1);
end

%% Solve for Steady state
disp('Solving Steady State by EGM')
tic
% Set parameters
defineSS_pars

mainskript_steadystate

VARk_SS=meshes.m(:).^2'*joint_distr(:)-(meshes.m(:)'*joint_distr(:))^2;
VARh_SS=meshes.h(:).^2'*joint_distr(:)-(meshes.h(:)'*joint_distr(:))^2;
COV_SS=(meshes.h(:).*meshes.m(:))'*joint_distr(:)-(meshes.h(:)'*joint_distr(:))*(meshes.m(:)'*joint_distr(:));
CORR_SS=COV_SS/(VARk_SS^0.5*VARh_SS^0.5); % Definition of corr

    
toc
%%
disp('Reiter method with state space reduction')
cd('One Asset HANC Capital (JPEG)')


tic
mainskript_2ndorder
toc




%% Produce IRFs and simulate
  % 2nd order IRF
  % eta0 : indicator of shock cell
  % x0 : initial condition of state variables
  % e: T (simulation periods) by ne, zeros for specific shock series

%   par.sigmaS=.07;
  sigma = 10*par.sigmaS;
  
  x0_simul=zeros(mpar.numstates,1);
  x0_simul(end)=sigma;
   
  e_simul=0*randn(mpar.maxlag,1);

    
  [Y_control_simul_2nd,X_state_simul_2nd] = simu_2nd(gx, hx, gxx, hxx, gss, hss, eta0, par.sigmaS, x0_simul, e_simul);
 
   
  X_state_simul_2nd=X_state_simul_2nd';
  Y_control_simul_2nd=Y_control_simul_2nd';
  IRF_state_sparse_simul_2nd=[X_state_simul_2nd;Y_control_simul_2nd];
  
   e_simul_ = [0;e_simul(1:end-1)];
  [Y_control_simul_1st,X_state_simul_1st] = simu_1st(gx, hx, eta0, mpar.maxlag ,x0_simul, e_simul_/par.sigmaS^(-1));
    
  X_state_simul_1st=X_state_simul_1st';
  Y_control_simul_1st=Y_control_simul_1st';
  IRF_state_sparse_simul_1st=[X_state_simul_1st;Y_control_simul_1st];
 
  cd ..

%% Plotting 
% First order
IRF_distr=Gamma_state*IRF_state_sparse_simul_1st(1:mpar.numstates-1,1:mpar.maxlag);

IRF_H=100*grid.h(1:end)*IRF_distr(mpar.nm+(1:mpar.nh),2:end)/par.H;
K=grid.m*IRF_distr((1:mpar.nm),1:end)+grid.K;
I=(K(2:end)-(1-par.delta)*K(1:end-1));
IRF_I=100*(I/(par.delta*grid.K)-1);
IRF_K=100*(K(1:end)./grid.K-1);
IRF_S=100*IRF_state_sparse_simul_1st(mpar.numstates-os+1,1:end-1);
Y=Output*(1+IRF_state_sparse_simul_1st(end-oc+2,1:end-1));
IRF_C=100*((Y-I)./(targets.Y-par.delta*grid.K)-1);
IRF_Y=100*IRF_state_sparse_simul_1st(end-oc+2,1:end-1);
IRF_Q=100*IRF_state_sparse_simul_1st(end-oc+1,1:end-1);
IRF_N=100*IRF_state_sparse_simul_1st(end-oc+5,1:end-1);
IRF_R=100*IRF_state_sparse_simul_1st(end-oc+4,1:end-1);

% Second order
IRF_state_sparse_simul_2nd(:,end)=[];

IRF_distr=Gamma_state*IRF_state_sparse_simul_2nd(1:mpar.numstates-1,1:mpar.maxlag);

IRF_H=100*grid.h(1:end)*IRF_distr(mpar.nm+(1:mpar.nh),2:end)/par.H;
K=grid.m*IRF_distr((1:mpar.nm),1:end)+grid.K;
I=(K(2:end)-(1-par.delta)*K(1:end-1));
IRF_I_2nd=100*(I/(par.delta*grid.K)-1);
IRF_K_2nd=100*(K(1:end)./grid.K-1);
IRF_S_2nd=100*IRF_state_sparse_simul_2nd(mpar.numstates-os+1,1:end-1);
Y=Output*(1+IRF_state_sparse_simul_2nd(end-oc+2,1:end-1));
IRF_C_2nd=100*((Y-I)./(targets.Y-par.delta*grid.K)-1);
IRF_Y_2nd=100*IRF_state_sparse_simul_2nd(end-oc+2,1:end-1);
IRF_Q_2nd=100*IRF_state_sparse_simul_2nd(end-oc+1,1:end-1);
IRF_N_2nd=100*IRF_state_sparse_simul_2nd(end-oc+5,1:end-1);
IRF_R_2nd=100*IRF_state_sparse_simul_2nd(end-oc+4,1:end-1);

%%
figurename=['IRF_K_ComparisonOrder_' casenamemaster  '_' aggrshock];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 1500 1000])
plot(1:mpar.maxlag,IRF_K,'b-','LineWidth',3.5)
hold on
plot(1:mpar.maxlag,IRF_K_2nd,'r--','LineWidth',3.5)
% ylim([0 0.25])
legend('1. Order Perturbation','2. Order Perturbation','Location','SouthEast')
ylabel('Percent','Interpreter','none','FontName','arial','FontSize',40)
xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',40)
set(gca, 'FontName','arial','FontSize',40);
set(gca,'Position',[0.125 0.15 0.85 0.825])
printpdf(gcf,['../results/figure6/' figurename])

%%
[Ey,Ex] = unconditional_mean(gx, hx, gxx, hxx, gss, hss, eta0, sigma);

EJD=Gamma_state*Ex(1:mpar.numstates-1);

EKD=sum(joint_distr,2)+EJD(1:mpar.nm);

figurename=['IRF_DISTR_ComparisonOrder_' casenamemaster  '_' aggrshock];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 1500 1000])
plot(grid.m./targets.K,sum(joint_distr,2),'b-','LineWidth',3.5)
hold on
plot(grid.m./targets.K,EKD,'r--','LineWidth',3.5)
% ylim([0 0.25])
legend('1. Order Perturbation','2. Order Perturbation','Location','NorthEast')
ylabel('Percent','Interpreter','none','FontName','arial','FontSize',40)
xlabel('Capital','Interpreter','none','FontName','arial','FontSize',40)
set(gca, 'FontName','arial','FontSize',40);
set(gca,'Position',[0.15 0.15 0.8 0.825])
printpdf(gcf,['../results/figure6/' figurename])

%%
close all