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
casenamemaster='SS_BASELINE_HANC_COP';
addpath(genpath('shared_functions'))
oldpath = path;


p = gcp('nocreate');
pool='local'; % select machine

if isempty(p)
    c = parcluster;
    parpool(pool,5)
    p = gcp('nocreate');
end

%% Select aggregate shock
aggrshock           = 'TFP';

par.rhoS            = 0.75;     % Persistence 
par.sigmaS          = 0.007;    % STD 
mpar.ns             = 3;        % number of TFP states
mpar.T(1,1)         = 1000;     % number of periods
mpar.T(2,1)         = 500;      % burn-in phase
mpar.maxlag         = 40;       % number of periods IRF

mpar.nK             = 3;        % number of grid points for KS

% Select aggregate shock
% aggrshock           = 'Uncertainty';
% par.rhoS            = 0.84;    % Persistence of variance
% par.sigmaS          = 0.54;    % STD of variance shocks


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
% defineSS_pars_UNC

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
mainskript
toc

cd ..
%% Produce IRFs and simulate
x0=zeros(mpar.numstates,1);
x0(end)=par.sigmaS;

MX=[eye(length(x0));gx];
IRF_state_sparse=[];
x=x0;

for t=1:mpar.maxlag
    IRF_state_sparse(:,t)=(MX*x)';
    x=hx*x;
end

IRF_distr=Gamma_state*IRF_state_sparse(1:mpar.numstates-1,1:mpar.maxlag);

IRF_H=100*grid.h(1:end)*IRF_distr(mpar.nm+(1:mpar.nh),2:end)/par.H;
K=grid.m*IRF_distr((1:mpar.nm),1:end)+grid.K;
I=(K(2:end)-(1-par.delta)*K(1:end-1));
IRF_I=100*(I/(par.delta*grid.K)-1);
IRF_K=100*(K(1:end)./grid.K-1);
IRF_S=100*IRF_state_sparse(mpar.numstates-os+1,1:end-1);
Y=Output*(1+IRF_state_sparse(end-oc+2,1:end-1));
IRF_C=100*((Y-I)./(targets.Y-par.delta*grid.K)-1);
IRF_Y=100*IRF_state_sparse(end-oc+2,1:end-1);
IRF_Q=100*IRF_state_sparse(end-oc+1,1:end-1);
IRF_N=100*IRF_state_sparse(end-oc+5,1:end-1);
IRF_R=100*IRF_state_sparse(end-oc+4,1:end-1);

save(['IRF_' casename],'IRF_K','IRF_I','IRF_S','IRF_Y','IRF_C','IRF_Q','IRF_N','IRF_R','IRF_H')

xlinear=x0;
IRF_state_sparse=[];
for t=1:mpar.T(1,1)
    xlinear(end)= log(grid.s(smaster(t)));
    IRF_state_sparse(:,t)=(MX*xlinear)';
    xlinear=hx*xlinear;
end

IRF_distr=Gamma_state*IRF_state_sparse(1:mpar.numstates-1,1:mpar.T(1,1));

Haux=(grid.h*IRF_distr(mpar.nm+(1:mpar.nh),1:end)+par.H)';
Kaux=(grid.m*IRF_distr((1:mpar.nm),1:end)+grid.K)';

Kmat(:,1)=Kaux;
Imat(:,1)=(Kaux(2:end)-(1-par.delta)*Kaux(1:end-1));



JDREDT=zeros(mpar.nm*mpar.nh,mpar.T(1,1));
JDREDT(:,1)=joint_distr(:);
for t=1:mpar.T(1,1)-1
    marginal_mminus =  Xss(1:mpar.nm)' + IRF_distr(1:mpar.nm,t)';
    marginal_hminus =  Xss(mpar.nm+(1:mpar.nh))' + IRF_distr(mpar.nm+(1:mpar.nh),t)';
    % Calculate joint distributions
    cumdist = zeros(mpar.nm+1,mpar.nh+1);
    cumdist(2:end,2:end) = Copula({cumsum(marginal_mminus),cumsum(marginal_hminus)});
    JDCOP = diff(diff(cumdist,1,1),1,2);
    JDREDT(:,t+1)=JDCOP(:)./sum(JDCOP(:));
        
    MEANk_RED(t)=meshes.m(:)'*JDCOP(:);
    VARk_RED(t)=meshes.m(:).^2'*JDCOP(:)-(meshes.m(:)'*JDCOP(:))^2;
    VARh_RED(t)=meshes.h(:).^2'*JDCOP(:)-(meshes.h(:)'*JDCOP(:))^2;
    COV_RED(t)=(meshes.h(:).*meshes.m(:))'*JDCOP(:)-(meshes.h(:)'*JDCOP(:))*(meshes.m(:)'*JDCOP(:));

    CORR_RED(t)=COV_RED(t)/(VARk_RED(t)^0.5*VARh_RED(t)^0.5); % Definition of corr
  
end

%%
disp('Reiter method with state space reduction plus perturbation of Copula')
cd('One Asset HANC Capital (JPEG-COPULA)')


tic
mainskript
toc

cd ..


%% Produce IRFs and simulate
x0=zeros(mpar.numstates,1);
x0(end)=par.sigmaS;

MX=[eye(length(x0));gx];
IRF_state_sparse=[];
x=x0;

for t=1:mpar.maxlag
    IRF_state_sparse(:,t)=(MX*x)';
    x=hx*x;
end

IRF_distr=Gamma_state*IRF_state_sparse(1:mpar.nm+mpar.nh-2,1:mpar.maxlag);

IRF_H=100*grid.h(1:end)*IRF_distr(mpar.nm+(1:mpar.nh),2:end)/par.H;
K=grid.m*IRF_distr((1:mpar.nm),1:end)+grid.K;
I=(K(2:end)-(1-par.delta)*K(1:end-1));
IRF_I=100*(I/(par.delta*grid.K)-1);
IRF_K=100*(K(1:end)./grid.K-1);
IRF_S=100*IRF_state_sparse(mpar.numstates-os+1,1:end-1);
Y=Output*(1+IRF_state_sparse(end-oc+2,1:end-1));
IRF_C=100*((Y-I)./(targets.Y-par.delta*grid.K)-1);
IRF_Y=100*IRF_state_sparse(end-oc+2,1:end-1);
IRF_Q=100*IRF_state_sparse(end-oc+1,1:end-1);
IRF_N=100*IRF_state_sparse(end-oc+5,1:end-1);
IRF_R=100*IRF_state_sparse(end-oc+4,1:end-1);

%% plot_IRFs

save(['IRF_' casename],'IRF_K','IRF_I','IRF_S','IRF_Y','IRF_C','IRF_Q','IRF_N','IRF_R','IRF_H')

xlinear=x0;
IRF_state_sparse=[];
for t=1:mpar.T(1,1)
    xlinear(end)= log(grid.s(smaster(t)));
    IRF_state_sparse(:,t)=(MX*xlinear)';
    xlinear=hx*xlinear;
end

IRF_distr=Gamma_state*IRF_state_sparse(1:mpar.nm+mpar.nh-2,1:mpar.T(1,1));

Haux=(grid.h*IRF_distr(mpar.nm+(1:mpar.nh),1:end)+par.H)';
Kaux=(grid.m*IRF_distr((1:mpar.nm),1:end)+grid.K)';

Kmat(:,2)=Kaux;
Imat(:,2)=(Kaux(2:end)-(1-par.delta)*Kaux(1:end-1));

COPind= mpar.nm+mpar.nh-2+(1:length(indexCOPdct));
NN   = mpar.nm*mpar.nh; % Number of points in the full grid
NN2   = (mpar.nm-1)*(mpar.nh-1); % Number of points in the full grid

JDCOPT=zeros(mpar.nm*mpar.nh,mpar.T(1,1));
JDCOPT(:,1)=joint_distr(:);
for t=1:mpar.T(1,1)-1
    marginal_mminus =  Xss(1:mpar.nm)' + IRF_distr(1:mpar.nm,t)';
    marginal_hminus =  Xss(mpar.nm+(1:mpar.nh))' + IRF_distr(mpar.nm+(1:mpar.nh),t)';
    % Calculate joint distributions
    
    % Controls copula
    XX             = zeros(NN2,1);
    XX(indexCOPdct) = IRF_state_sparse(COPind);
    cum_aux = reshape(XX,[mpar.nm-1, mpar.nh-1]);
    cum_aux = myidct(cum_aux,1,DC3); % do dct-transformation
    cum_aux = myidct(cum_aux,2,DC4); % do dct-transformation
    cum_aux2= zeros(mpar.nm,mpar.nh);
    cum_aux2(1:end-1,1:end-1)=cum_aux;
    cum_aux2(end,1:end-1) = -sum(cum_aux,1); 
    cum_aux2(:,end) = -sum(cum_aux2,2); 
    cum_aux = (cum_aux2(:)+Xss((mpar.nm+mpar.nh)+(1:NN)));
    cum_aux = max(reshape(cum_aux,[mpar.nm, mpar.nh]),1e-16);
    cum_aux = (cum_aux)./sum(cum_aux(:));
    cum_aux = cumsum(cumsum(cum_aux,1),2);
    
    cum_zero=zeros(mpar.nm+1,mpar.nh+1);
    cum_zero(2:end,2:end)=cum_aux;
    Copula = griddedInterpolant({cum_zero(:,end),cum_zero(end,:)},cum_zero,'spline');

    cumdist = zeros(mpar.nm+1,mpar.nh+1);
    cumdist(2:end,2:end) = Copula({cumsum(marginal_mminus),cumsum(marginal_hminus)});
    JDCOP = diff(diff(cumdist,1,1),1,2);
    JDCOPT(:,t+1)=JDCOP(:);
    
    MEANk_CP(t)=meshes.m(:)'*JDCOP(:);
    VARk_CP(t)=meshes.m(:).^2'*JDCOP(:)-(meshes.m(:)'*JDCOP(:))^2;
    VARk_CP2(t)=grid.m(:).^2'*marginal_mminus(:)-(grid.m(:)'*marginal_mminus(:))^2;
    VARh_CP(t)=meshes.h(:).^2'*JDCOP(:)-(meshes.h(:)'*JDCOP(:))^2;
    COV_CP(t)=(meshes.h(:).*meshes.m(:))'*JDCOP(:)-(meshes.h(:)'*JDCOP(:))*(meshes.m(:)'*JDCOP(:));

    CORR_CP(t)=COV_CP(t)/(VARk_CP(t)^0.5*VARh_CP(t)^0.5); % Definition of corr
  
end

%%
disp('Reiter method without reduction')
cd('One Asset HANC Capital (FULL)')

tic
mainskript
toc

cd ..
%% Produce IRFs and simulate
x0=zeros(mpar.numstates,1);
x0(end)=par.sigmaS;

MX=[eye(length(x0));gx];
IRF_state_sparse=[];
x=x0;

for t=1:mpar.maxlag
    IRF_state_sparse(:,t)=(MX*x)';
    x=hx*x;
end


IRF_distr=Gamma_state*IRF_state_sparse(1:mpar.numstates-os,1:mpar.maxlag);

IRF_distr=reshape(IRF_distr,[mpar.nm mpar.nh mpar.maxlag]);
IRF_distr_m=squeeze(sum(IRF_distr,2));
IRF_distr_h=squeeze(sum(IRF_distr,1));
IRF_H=100*grid.h(1:end)*IRF_distr_h(1:end,2:end)/par.H;
K=grid.m*IRF_distr_m(:,1:end)+grid.K;
I=(K(2:end)-(1-par.delta)*K(1:end-1));
IRF_I=100*(I/(par.delta*grid.K)-1);
IRF_K=100*(K(1:end)./grid.K-1);
IRF_S=100*IRF_state_sparse(mpar.numstates-os+1,1:end-1);
Y=Output*(1+IRF_state_sparse(end-oc+2,1:end-1));
IRF_C=100*((Y-I)./(targets.Y-par.delta*grid.K)-1);
IRF_Y=100*IRF_state_sparse(end-oc+2,1:end-1);
IRF_Q=100*IRF_state_sparse(end-oc+1,1:end-1);
IRF_N=100*IRF_state_sparse(end-oc+5,1:end-1);
IRF_R=100*IRF_state_sparse(end-oc+4,1:end-1);

% plot_IRFs

save(['IRF_' casename],'IRF_K','IRF_I','IRF_S','IRF_Y','IRF_C','IRF_Q','IRF_N','IRF_R','IRF_H')

xlinear=x0;
MX=[eye(length(x0));gx];
IRF_state_sparse=[];
for t=1:mpar.T(1,1)
    xlinear(end)= log(grid.s(smaster(t)));
    IRF_state_sparse(:,t)=(MX*xlinear)';
    xlinear=hx*xlinear;
end

IRF_distr=Gamma_state*IRF_state_sparse(1:mpar.numstates-os,1:mpar.T(1,1));

IRF_distr=reshape(IRF_distr,[mpar.nm mpar.nh mpar.T(1,1)]);
IRF_distr_m=squeeze(sum(IRF_distr,2));
IRF_distr_h=squeeze(sum(IRF_distr,1));
Kaux=grid.m*IRF_distr_m(:,1:end)+grid.K;
Kmat(:,3)=Kaux;
Imat(:,3)=(Kaux(2:end)-(1-par.delta)*Kaux(1:end-1));


JDTR=zeros(mpar.nm*mpar.nh,mpar.T(1,1));
JDTR(:,1)=joint_distr(:);
for t=1:mpar.T(1,1)-1
    JDaux =  joint_distr + IRF_distr(:,:,t);
    
    JDTR(:,t+1)=JDaux(:);
    
    MEANk_R(t)=meshes.m(:)'*JDaux(:);

    VARk_R(t)=meshes.m(:).^2'*JDaux(:)-(meshes.m(:)'*JDaux(:))^2;
    VARh_R(t)=meshes.h(:).^2'*JDaux(:)-(meshes.h(:)'*JDaux(:))^2;
    COV_R(t)=(meshes.h(:).*meshes.m(:))'*JDaux(:)-(meshes.h(:)'*JDaux(:))*(meshes.m(:)'*JDaux(:));

    CORR_R(t)=COV_R(t)/(VARk_R(t)^0.5*VARh_R(t)^0.5); % Definition of corr
  
end


%%

for t=1:mpar.T(1,1)
    AA = sum(JDTR(:,t).*log(2*JDTR(:,t)./(JDCOPT(:,t)+JDTR(:,t))));
    BB = sum(JDCOPT(:,t).*log(2*JDCOPT(:,t)./(JDCOPT(:,t)+JDTR(:,t))));
    RD_COP(t)=(AA+BB)/2;
end

for t=1:mpar.T(1,1)
    AA = sum(JDTR(:,t).*log(2*JDTR(:,t)./(JDREDT(:,t)+JDTR(:,t))));
    BB = sum(JDREDT(:,t).*log(2*JDREDT(:,t)./(JDREDT(:,t)+JDTR(:,t))));
    RD_RED(t)=(AA+BB)/2;
end

for t=1:mpar.T(1,1)
    aux =(JDTR(:,t)+joint_distr(:));
    AA = sum(joint_distr(:).*log(2*joint_distr(:)./(aux)));
    BB = sum(JDTR(:,t).*log(2*JDTR(:,t)./(aux)));
    RD_SS(t)=(AA+BB)/2;
end
%%
figurename=['DISTR_Comparison_' casenamemaster  '_' aggrshock];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 1500 1000])
plot(mpar.T(2,1)+1:mpar.T(1,1),real(RD_RED(mpar.T(2,1)+1:end).^0.5),'r-','LineWidth',3.5)
hold on
plot(mpar.T(2,1)+1:mpar.T(1,1),real(RD_COP(mpar.T(2,1)+1:end).^0.5),'b--','LineWidth',3.5)
% ylim([0 1.5e-3])
ylabel('Jensen-Shannon Distance','Interpreter','none','FontName','arial','FontSize',40)
xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',40)
set(gca, 'FontName','arial','FontSize',40);
legend('Reiter-Reduction','Reiter-Reduction w Copula Perturbation')
set(gca,'Position',[0.15 0.15 0.8 0.75])
printpdf(gcf,['../results/figure5/' figurename])

figurename=['DISTR_ComparisonALL_' casenamemaster  '_' aggrshock];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 1500 1000])
plot(mpar.T(2,1)+1:mpar.T(1,1),real(RD_RED(mpar.T(2,1)+1:end).^0.5),'b-','LineWidth',3.5)
hold on
plot(mpar.T(2,1)+1:mpar.T(1,1),real(RD_COP(mpar.T(2,1)+1:end).^0.5),'r--','LineWidth',3.5)
ylim([0 2e-2])
plot(mpar.T(2,1)+1:mpar.T(1,1),real(RD_SS(mpar.T(2,1)+1:end).^0.5),'k:','LineWidth',3.5)

ylabel('Jensen-Shannon Distance','Interpreter','none','FontName','arial','FontSize',40)
xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',40)
set(gca, 'FontName','arial','FontSize',40);
legend('Reiter-Red. (vs Reiter-Full)','Reiter-Red. w Copula Pert. (vs Reiter-Full)', 'Reiter-Full (vs SS)')
set(gca,'Position',[0.175 0.15 0.775 0.8])
printpdf(gcf,['../results/figure5/' figurename])

figurename=['DISTR_ComparisonSS_' casenamemaster  '_' aggrshock];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 1500 1000])
plot(mpar.T(2,1)+1:mpar.T(1,1),real(RD_RED(mpar.T(2,1)+1:end).^0.5),'-','LineWidth',3.5)
hold on
plot(mpar.T(2,1)+1:mpar.T(1,1),real(RD_SS(mpar.T(2,1)+1:end).^0.5),'--','LineWidth',3.5)
% ylim([0 1.5e-3])
ylabel('Jensen-Shannon Distance','Interpreter','none','FontName','arial','FontSize',40)
xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',40)
set(gca, 'FontName','arial','FontSize',40);
legend('Reiter-Reduction vs Reiter-Full','Reiter-Reduction vs SS')
set(gca,'Position',[0.15 0.15 0.8 0.75])
printpdf(gcf,['../results/figure5/' figurename])

ErrorTabletoDISTNorm(1,:)=mean((real([[RD_RED(mpar.T(2,1)+1:end); RD_COP(mpar.T(2,1)+1:end)].^0.5])),2)';
ErrorTabletoDISTNorm(2,:)=real([max(RD_RED(mpar.T(2,1)+1:end).^0.5); max(RD_COP(mpar.T(2,1)+1:end).^0.5)]');

filename=['../results/AvgAbsErrorDist_' casenamemaster '_' aggrshock '.tex'];
FID = fopen(filename, 'w');
fprintf(FID, '   %8.4f \n', ErrorTabletoDISTNorm(1,1:1));
fclose(FID);
filename=['../results/MaxAbsErrorDist_' casenamemaster '_' aggrshock '.tex'];
FID = fopen(filename, 'w');
fprintf(FID, '   %8.4f \n', ErrorTabletoDISTNorm(2,1:1));
fclose(FID);

figurename=['K_Comparison_' casenamemaster  '_' aggrshock];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 1500 1000])
plot(mpar.T(2,1)+1:mpar.T(1,1),Kmat(mpar.T(2,1)+1:end,1),'b-','LineWidth',3.5)
hold on
plot(mpar.T(2,1)+1:mpar.T(1,1),Kmat(mpar.T(2,1)+1:end,2),'r--','LineWidth',3.5)
plot(mpar.T(2,1)+1:mpar.T(1,1),Kmat(mpar.T(2,1)+1:end,3),'k:','LineWidth',3.5)
ylim([34 36])

% ylim([31.2 31.6])
ylabel('Capital','Interpreter','none','FontName','arial','FontSize',40)
xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',40)
set(gca, 'FontName','arial','FontSize',40);
legend('Reiter-Red.','Reiter-Red. w Copula Pert.', 'Reiter-Full','Location','NorthWest')
set(gca,'Position',[0.175 0.15 0.775 0.8])
printpdf(gcf,['../results/figure5/' figurename])

%%
close all
