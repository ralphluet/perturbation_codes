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
casenamemaster='SS_BASELINE_HANC';
addpath(genpath('shared_functions'))
oldpath = path;


p = gcp('nocreate');
pool='local'; % select machine

if isempty(p)
    c = parcluster;
%     parpool(pool,c.NumWorkers)
    parpool(pool,2)

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
mainskript
toc

%% Simulations
% DenHaan error
JDminus=joint_distr;
JDT=zeros(mpar.nm*mpar.nh,mpar.T(1,1));
JDT(:,1)=JDminus(:);
Krealized_RC=targets.K;
Hrealized_RC=par.H;
x=zeros(mpar.numstates,1);

PIGamma= pinv(Gamma_state);

for t=1:mpar.T(1,1)-1
    ControlNOW=(gx*x)';
    xFOR=hx*x;
    ControlFOR=(gx*xFOR)';
    
    [JDminus]  = denHaan_Reiter(JDminus,log(grid.s(smaster(t))),...
        ControlFOR',  ControlNOW',...
        Yss,indexMUdct,par,mpar,grid,meshes,P_H,aggrshock,oc);

    Krealized_RC(t+1)=meshes.m(:)'*JDminus(:);
    Hrealized_RC(t+1)=meshes.h(:)'*JDminus(:);

    VARk_DH(t)=meshes.m(:).^2'*JDminus(:)-(meshes.m(:)'*JDminus(:))^2;
    VARh_DH(t)=meshes.h(:).^2'*JDminus(:)-(meshes.h(:)'*JDminus(:))^2;
    COV_DH(t)=(meshes.h(:).*meshes.m(:))'*JDminus(:)-(meshes.h(:)'*JDminus(:))*(meshes.m(:)'*JDminus(:));

    CORR_DH(t)=COV_DH(t)/(VARk_DH(t)^0.5*VARh_DH(t)^0.5); % Definition of corr
    
    Hrealized_RC(t+1)=meshes.h(:)'*JDminus(:);
    % liquid assets
    aux_m = squeeze(sum(JDminus,2));
    % human capital
    aux_h = squeeze(sum(JDminus,1))';
    Xaux=([aux_m(1:end);aux_h(1:end)]-Xss(1:(mpar.nm+mpar.nh)));
    
    JDT(:,t+1)= JDminus(:);
    x(1:end-1)=PIGamma*Xaux;
    x(end)=log(grid.s(smaster(t+1)));    
end
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

% plot_IRFs

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

DenHaanError(:,1)=100*log(Kaux./Krealized_RC');



JDCOPT=zeros(mpar.nm*mpar.nh,mpar.T(1,1));
JDCOPT(:,1)=joint_distr(:);
for t=1:mpar.T(1,1)-1
    marginal_mminus =  Xss(1:mpar.nm)' + IRF_distr(1:mpar.nm,t)';
    marginal_hminus =  Xss(mpar.nm+(1:mpar.nh))' + IRF_distr(mpar.nm+(1:mpar.nh),t)';
    % Calculate joint distributions
    cumdist = zeros(mpar.nm+1,mpar.nh+1);
    cumdist(2:end,2:end) = Copula({cumsum(marginal_mminus),cumsum(marginal_hminus)});
    JDCOP = diff(diff(cumdist,1,1),1,2);
    JDCOPT(:,t+1)=JDCOP(:);
    
    VARk_CP(t)=meshes.m(:).^2'*JDCOP(:)-(meshes.m(:)'*JDCOP(:))^2;
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

%% Simulations
% DenHaan error
JDminus=joint_distr;
Krealized_RF=targets.K;
Hrealized_RF=par.H;
x=zeros(mpar.numstates,1);

PIGamma= pinv(Gamma_state);

for t=1:mpar.T(1,1)-1
    ControlNOW=(gx*x)';
    xFOR=hx*x;
    ControlFOR=(gx*xFOR)';
    
    [JDminus]  = denHaan_Reiter(JDminus,log(grid.s(smaster(t))),...
        ControlFOR',  ControlNOW',...
        Yss,indexMUdct,par,mpar,grid,meshes,P_H,aggrshock,oc);

    Krealized_RF(t+1)=meshes.m(:)'*JDminus(:);
    Hrealized_RF(t+1)=meshes.h(:)'*JDminus(:);

    Xaux=([JDminus(:)]-Xss(1:(mpar.nm*mpar.nh)));
    
    x(1:end-1)=PIGamma*Xaux;
    x(end)=log(grid.s(smaster(t+1)));    
end
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
Kmat(:,2)=Kaux;
Imat(:,2)=(Kaux(2:end)-(1-par.delta)*Kaux(1:end-1));

DenHaanError(:,2)=100*log(Kaux./Krealized_RF);


JDTR=zeros(mpar.nm*mpar.nh,mpar.T(1,1));
JDTR(:,1)=joint_distr(:);
for t=1:mpar.T(1,1)-1
    JDaux =  joint_distr + IRF_distr(:,:,t);
    
    JDTR(:,t+1)=JDaux(:);
    
    VARk_R(t)=meshes.m(:).^2'*JDaux(:)-(meshes.m(:)'*JDaux(:))^2;
    VARh_R(t)=meshes.h(:).^2'*JDaux(:)-(meshes.h(:)'*JDaux(:))^2;
    COV_R(t)=(meshes.h(:).*meshes.m(:))'*JDaux(:)-(meshes.h(:)'*JDaux(:))*(meshes.m(:)'*JDaux(:));

    CORR_R(t)=COV_R(t)/(VARk_R(t)^0.5*VARh_R(t)^0.5); % Definition of corr
  
end


for t=1:mpar.T(1,1)
    AA = sum(JDTR(:,t).*log(2*JDTR(:,t)./(JDCOPT(:,t)+JDTR(:,t))));
    BB = sum(JDCOPT(:,t).*log(2*JDCOPT(:,t)./(JDCOPT(:,t)+JDTR(:,t))));
    RD(t)=(AA+BB)/2;
end



figurename=['DISTR_Comparison_' casenamemaster  '_' aggrshock];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 1500 1000])
plot(1:mpar.T(1,1)-1,RD(2:end).^0.5,'LineWidth',3.5)
ylabel('Jensen-Shannon Distance','Interpreter','none','FontName','arial','FontSize',40)
xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',40)
set(gca, 'FontName','arial','FontSize',40);
set(gca,'Position',[0.15 0.15 0.8 0.75])
printpdf(gcf,['../results/' figurename])

ErrorTabletoDISTNorm(1,:)=mean(abs([RD.^0.5]));
ErrorTabletoDISTNorm(2,:)=max(abs([RD.^0.5]));

filename=['../results/AvgAbsErrorDist_' casenamemaster '_' aggrshock '.tex'];
FID = fopen(filename, 'w');
fprintf(FID, '   %8.4f \n', ErrorTabletoDISTNorm(1,1:1));
fclose(FID);
filename=['../results/MaxAbsErrorDist_' casenamemaster '_' aggrshock '.tex'];
FID = fopen(filename, 'w');
fprintf(FID, '   %8.4f \n', ErrorTabletoDISTNorm(2,1:1));
fclose(FID);

%%
disp('Krusell Smith method')
cd('One Asset HANC Capital (KS)')

tic
mainskript
toc
cd ..
addpath(genpath('One Asset HANC Capital (KS)/functions'))

[IRF_K_KS, IRF_S_KS]=KSsimulationIRF(k_star,P_H,joint_distr,grid,mesh,mpar,par);

[Kmat(:,3), ss_aux]=KSsimulation(k_star,P_H,PS_S,pr_s_master(:),joint_distr,grid,mesh,mpar);
Kaux=Kmat(:,3);
Imat(:,3)=(Kaux(2:end)-(1-par.delta)*Kaux(1:end-1));

Klinear_KS=zeros(mpar.T(1,1),1);
Klinear_KS(1)=targets.K;
for t=1:mpar.T(1,1)-1
    
    Klinear_KS(t+1)=exp(par.beta_K(1,smaster(t))   + par.beta_K(2,smaster(t)).*log(Klinear_KS(t)));
%             Klinear_KS(t+1)     = exp(par.beta_K(1)   + par.beta_K(2)*log(grid.s(smaster(t))) + par.beta_K(3)*log(Klinear_KS(t)));

    
end

DenHaanError(:,3)=100*log(Kaux./Klinear_KS);

mean(abs(DenHaanError),1)

path(oldpath)
%% Save simulation results
figurename=['K_Comparison_' casenamemaster  '_' aggrshock];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 1500 1000])
plot(mpar.T(2,1)+1:mpar.T(1,1),Kmat(mpar.T(2,1)+1:end,1),'b-','LineWidth',3.5)
hold on
plot(mpar.T(2,1)+1:mpar.T(1,1),Kmat(mpar.T(2,1)+1:end,2),'r--','LineWidth',3.5)
plot(mpar.T(2,1)+1:mpar.T(1,1),Kmat(mpar.T(2,1)+1:end,3),'k:','LineWidth',3.5)
legend('Reiter-Reduction','Reiter-Full','Krusell-Smith','Location','SouthEast')
ylabel('Capital','Interpreter','none','FontName','arial','FontSize',40)
xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',40)
set(gca, 'FontName','arial','FontSize',40);
set(gca,'Position',[0.15 0.15 0.8 0.825])
printpdf(gcf,['../results/figure1/' figurename])

figurename=['K_ComparisonShort_' casenamemaster  '_' aggrshock];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 1500 1000])
plot(mpar.T(1,1)-100:mpar.T(1,1),Kmat(mpar.T(1,1)-100:end,1),'b-','LineWidth',3.5)
hold on
plot(mpar.T(1,1)-100:mpar.T(1,1),Kmat(mpar.T(1,1)-100:end,2),'r--','LineWidth',3.5)
plot(mpar.T(1,1)-100:mpar.T(1,1),Kmat(mpar.T(1,1)-100:end,3),'k:','LineWidth',3.5)
legend('Reiter-Reduction','Reiter-Full','Krusell-Smith','Location','SouthEast')
ylabel('Capital','Interpreter','none','FontName','arial','FontSize',40)
xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',40)
set(gca, 'FontName','arial','FontSize',40);
set(gca,'Position',[0.15 0.15 0.8 0.825])
printpdf(gcf,['../results/figure1/' figurename])

Kmat_normalized=Kmat./repmat(mean(Kmat(mpar.T(2,1):end,:),1),[mpar.T(1,1) 1]);
RCtoKS=100*log(Kmat_normalized(mpar.T(2,1):end,1)./Kmat_normalized(mpar.T(2,1):end,3));
RFtoKS=100*log(Kmat_normalized(mpar.T(2,1):end,2)./Kmat_normalized(mpar.T(2,1):end,3));

ErrorTabletoKSNorm(1,:)=mean(abs([RCtoKS RFtoKS]));
ErrorTabletoKSNorm(2,:)=max(abs([RCtoKS RFtoKS]));


RCtoKS=100*log(Kmat(:,1)./Kmat(:,3));
RFtoKS=100*log(Kmat(:,2)./Kmat(:,3));

ErrorTabletoKS(1,:)=mean(abs([RCtoKS RFtoKS]));
ErrorTabletoKS(2,:)=max(abs([RCtoKS RFtoKS]));


filename=['../results/table1/AvgAbsErrortoKS_' casenamemaster '_' aggrshock '.tex'];
FID = fopen(filename, 'w');
fprintf(FID, '   %8.4f & %8.4f \n', ErrorTabletoKS(1,1:2));
fclose(FID);
filename=['../results/table1/MaxAbsErrortoKS_' casenamemaster '_' aggrshock '.tex'];
FID = fopen(filename, 'w');
fprintf(FID, '   %8.4f & %8.4f \n', ErrorTabletoKS(2,1:2));
fclose(FID);

RCtoRF=100*log(Kmat(:,1)./Kmat(:,2));
ErrorTabletoRF(1,:)=mean(abs([RCtoRF]));
ErrorTabletoRF(2,:)=max(abs([RCtoRF]));

filename=['../results/table1/AvgAbsErrortoRF_' casenamemaster '_' aggrshock '.tex'];
FID = fopen(filename, 'w');
fprintf(FID, '   %8.4f \n', ErrorTabletoRF(1,:));
fclose(FID);
filename=['../results/table1/MaxAbsErrortoRF_' casenamemaster '_' aggrshock '.tex'];
FID = fopen(filename, 'w');
fprintf(FID, '   %8.4f \n', ErrorTabletoRF(2,:));
fclose(FID);

filename=['../results/table2/AvgAbsErrordenHaan_' casenamemaster '_' aggrshock '.tex'];
FID = fopen(filename, 'w');
fprintf(FID, '   %8.4f & %8.4f & %8.4f \n', mean(abs(DenHaanError),1));
fclose(FID);
filename=['../results/table2/MaxAbsErrordenHaan_' casenamemaster '_' aggrshock '.tex'];
FID = fopen(filename, 'w');
fprintf(FID, '   %8.4f & %8.4f & %8.4f \n',max(abs(DenHaanError),[],1));
fclose(FID);


%% Plot policy functions and IRFs for various number of basis functions
disp('Compare approximation quality')
% Set parameters
defineSS_pars

mainskript_steadystate

aux = reshape(c_guess,[mpar.nm, mpar.nh]);
aux = mydct(mydct(aux,1),2); % do dct-transformation

THETA=aux(:);
[~,ind] = sort(abs(THETA(:)),'descend');

for keepnumber=[10 50 200];

compressionIndexes=ind(1:keepnumber);
theta = zeros(size(c_guess));
theta = theta(:);
theta(compressionIndexes) = THETA(compressionIndexes);
aux = reshape(theta,[mpar.nm, mpar.nh]);
c_reduced = myidct(myidct(aux,1),2); % do dct-transformation

figurename=['IRF_Policy_Comparison_' casenamemaster '_noBF_' num2str(keepnumber)];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 1500 1000])
plot(grid.m(1:end-40),c_guess(1:end-40,1),'r--','LineWidth',3.5)
hold on
plot(grid.m(1:end-40),c_guess(1:end-40,2),'b--','LineWidth',3.5)

plot(grid.m(1:end-40),c_reduced(1:end-40,1),'r','LineWidth',3.5)
plot(grid.m(1:end-40),c_reduced(1:end-40,2),'b','LineWidth',3.5)

legend('Full grid, low income','Full grid, high income','Approximation, low income','Approximation, high income','Location','SouthEast')
ylabel('Consumption','Interpreter','none','FontName','arial','FontSize',40)
xlabel('Capital','Interpreter','none','FontName','arial','FontSize',40)
set(gca, 'FontName','arial','FontSize',40);
set(gca,'Position',[0.125 0.15 0.85 0.825])
printpdf(gcf,['../results/figure2/' figurename])

% cd('One Asset HANC Capital (JPEG)')
addpath(genpath('One Asset HANC Capital (JPEG)'))

mainskript_plotting

% cd ..
path(oldpath)

x0=zeros(mpar.numstates,1);
x0(end)=0.2;%par.sigmaS;

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
IRF_K=100*(K(1:end-1)./grid.K-1);
IRF_S=100*IRF_state_sparse(mpar.numstates-os+1,1:end-1);
Y=Output*(1+IRF_state_sparse(end-oc+2,1:end-1));
IRF_C=100*((Y-I)./(targets.Y-par.delta*grid.K)-1);
IRF_Y=100*IRF_state_sparse(end-oc+2,1:end-1);
IRF_Q=100*IRF_state_sparse(end-oc+1,1:end-1);
IRF_N=100*IRF_state_sparse(end-oc+5,1:end-1);
IRF_R=100*IRF_state_sparse(end-oc+4,1:end-1);

figurename=['IRF_IRF_Comparison_' casenamemaster '_noBF_' num2str(keepnumber)];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 1500 1000])
plot(1:mpar.maxlag-1,IRF_I,'r--','LineWidth',3.5)
hold on
plot(1:mpar.maxlag-1,IRF_K,'b-','LineWidth',3.5)
plot(1:mpar.maxlag-1,IRF_Y,'ko','LineWidth',3.5)
plot(1:mpar.maxlag-1,IRF_C,'g-','LineWidth',3.5)
legend('Investment','Capital','Output','Consumption','Location','NorthEast')
ylabel('Percent','Interpreter','none','FontName','arial','FontSize',40)
xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',40)
set(gca, 'FontName','arial','FontSize',40);
ylim([-10 60]);
set(gca,'Position',[0.125 0.15 0.85 0.825])
printpdf(gcf,['../results/figure4/' figurename])

    DC1  = mydctmx(mpar.nm);
    DC2  = mydctmx(mpar.nh);
% Controls
XX             = zeros(mpar.nm*mpar.nh,1);
XX(indexMUdct) = IRF_state_sparse(mpar.numstates+[1:keepnumber],1);
aux = reshape(XX,[mpar.nm, mpar.nh]);
c_redux=DC1'*aux*DC2;
c_redux = c_redux(:)+Yss(1:mpar.nm*mpar.nh);
c_redux = reshape(c_redux,[mpar.nm mpar.nh]);

figurename=['IRF_PolicyDev_Comparison_' casenamemaster '_noBF_' num2str(keepnumber)];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 1500 1000])
plot(grid.m(1:end-40),(c_guess(1:end-40,1)),'r--','LineWidth',3.5)
hold on
plot(grid.m(1:end-40),(c_guess(1:end-40,2)),'b--','LineWidth',3.5)

plot(grid.m(1:end-40),(c_redux(1:end-40,1)),'r','LineWidth',3.5)
plot(grid.m(1:end-40),(c_redux(1:end-40,2)),'b','LineWidth',3.5)

legend('Full grid, low income','Full grid, high income','low income, +20% TFP-Shock','high income, +20% TFP-Shock','Location','SouthEast')
ylabel('Consumption','Interpreter','none','FontName','arial','FontSize',40)
xlabel('Capital','Interpreter','none','FontName','arial','FontSize',40)
set(gca, 'FontName','arial','FontSize',40);
set(gca,'Position',[0.125 0.15 0.85 0.825])
printpdf(gcf,['../results/figure3/' figurename])

end

%% Plot policy functions and IRFs for various number of basis functions
disp('Compare approximation quality')
clear MAXERRORPOLY MEANERRORPOLY Compression KT

% Set parameters
defineSS_pars

mainskript_steadystate

aux = reshape(c_guess,[mpar.nm, mpar.nh]);
aux = mydct(mydct(aux,1),2); % do dct-transformation
% aux = mydct(aux,2); % do dct-transformation

THETA=aux(:);
[~,ind] = sort(abs(THETA(:)),'descend');
count=0;

for DEGREEN=[200 50 40 30 20];
count=count+1;

% cheby

% cd('One Asset HANC Capital (JPEG)')
addpath(genpath('One Asset HANC Capital (JPEG)'))

mainskript_plotting_cheby

% cd ..
path(oldpath)

x0=zeros(mpar.numstates,1);
x0(end)=0.2;%par.sigmaS;

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
IRF_K=100*(K(1:end-1)./grid.K-1);
IRF_S=100*IRF_state_sparse(mpar.numstates-os+1,1:end-1);
Y=Output*(1+IRF_state_sparse(end-oc+2,1:end-1));
IRF_C=100*((Y-I)./(targets.Y-par.delta*grid.K)-1);
IRF_Y=100*IRF_state_sparse(end-oc+2,1:end-1);
IRF_Q=100*IRF_state_sparse(end-oc+1,1:end-1);
IRF_N=100*IRF_state_sparse(end-oc+5,1:end-1);
IRF_R=100*IRF_state_sparse(end-oc+4,1:end-1);

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

KT(:,2,count)=Kaux;

    DC1  = mydctmx(mpar.nm);
    DC2  = mydctmx(mpar.nh);
% Controls
keepnumber=length(indexMUdct);
keepN(count)=keepnumber;

XX             = zeros(mpar.nm*mpar.nh,1);
XX(indexMUdct) = IRF_state_sparse(mpar.numstates+[1:keepnumber],1);
aux = reshape(XX,[mpar.nm, mpar.nh]);
c_redux=DC1'*aux*DC2;
c_redux = c_redux(:)+Yss(1:mpar.nm*mpar.nh);
c_redux = reshape(c_redux,[mpar.nm mpar.nh]);

% dct
% cd('One Asset HANC Capital (JPEG)')
addpath(genpath('One Asset HANC Capital (JPEG)'))

mainskript_plotting

% cd ..
path(oldpath)

compressionIndexes=ind(1:keepnumber);
theta = zeros(size(c_guess));
theta = theta(:);
theta(compressionIndexes) = THETA(compressionIndexes);
aux = reshape(theta,[mpar.nm, mpar.nh]);
c_reduced = myidct(myidct(aux,1),2); % do dct-transformation


x0=zeros(mpar.numstates,1);
x0(end)=0.2;%par.sigmaS;

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
IRF_K=100*(K(1:end-1)./grid.K-1);
IRF_S=100*IRF_state_sparse(mpar.numstates-os+1,1:end-1);
Y=Output*(1+IRF_state_sparse(end-oc+2,1:end-1));
IRF_C=100*((Y-I)./(targets.Y-par.delta*grid.K)-1);
IRF_Y=100*IRF_state_sparse(end-oc+2,1:end-1);
IRF_Q=100*IRF_state_sparse(end-oc+1,1:end-1);
IRF_N=100*IRF_state_sparse(end-oc+5,1:end-1);
IRF_R=100*IRF_state_sparse(end-oc+4,1:end-1);


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

KT(:,1,count)=Kaux;

    DC1  = mydctmx(mpar.nm);
    DC2  = mydctmx(mpar.nh);
% Controls
XX             = zeros(mpar.nm*mpar.nh,1);
XX(indexMUdct) = IRF_state_sparse(mpar.numstates+[1:keepnumber],1);
aux = reshape(XX,[mpar.nm, mpar.nh]);
c_redux=DC1'*aux*DC2;

c_redux = c_redux(:)+Yss(1:mpar.nm*mpar.nh);
c_redux = reshape(c_redux,[mpar.nm mpar.nh]);

Compression(1,count)=DEGREEN;
Compression(2,count)=keepnumber;


end
%%
for tt=1:count
    MAXERRORPOLY(1,tt)=max(abs(100*(log(KT(:,1,tt)./KT(:,1,1)))));
    MAXERRORPOLY(2,tt)=max(abs(100*(log(KT(:,2,tt)./KT(:,2,1)))));
    MEANERRORPOLY(1,tt)=mean(abs(100*(log(KT(:,1,tt)./KT(:,1,1)))));
    MEANERRORPOLY(2,tt)=mean(abs(100*(log(KT(:,2,tt)./KT(:,2,1)))));
end

MAXERRORPOLY=100*100*100*MAXERRORPOLY;
MEANERRORPOLY=100*100*100*MEANERRORPOLY;


filename=['../results/table4/MaxAbsErrorPOLYFULL_' casenamemaster '_' aggrshock '.tex'];
FID = fopen(filename, 'w');
fprintf(FID, '   %8.4f & %8.4f & %8.4f &  %8.4f & %8.4f  \n',MAXERRORPOLY(1,:));
% fprintf(FID, '   %8.4f & %8.4f & %8.4f &  %8.4f & %8.4f & %8.4f \n',MAXERRORPOLY(2,:));
fclose(FID);

filename=['../results/table4/MaxAbsErrorPOLYDCT_' casenamemaster '_' aggrshock '.tex'];
FID = fopen(filename, 'w');
% fprintf(FID, '   %8.4f & %8.4f & %8.4f &  %8.4f & %8.4f & %8.4f \n',MAXERRORPOLY(1,:));
fprintf(FID, '   %8.4f & %8.4f & %8.4f &  %8.4f & %8.4f  \n',MAXERRORPOLY(2,:));
fclose(FID);

filename=['../results/table4/MeanAbsErrorPOLYFULL_' casenamemaster '_' aggrshock '.tex'];
FID = fopen(filename, 'w');
fprintf(FID, '   %8.4f & %8.4f & %8.4f &  %8.4f & %8.4f  \n',MEANERRORPOLY(1,:));
% fprintf(FID, '   %8.4f & %8.4f & %8.4f &  %8.4f & %8.4f & %8.4f \n',MEANERRORPOLY(2,:));
fclose(FID);

filename=['../results/table4/MeanAbsErrorPOLYDCT_' casenamemaster '_' aggrshock '.tex'];
FID = fopen(filename, 'w');
% fprintf(FID, '   %8.4f & %8.4f & %8.4f &  %8.4f & %8.4f & %8.4f \n',MEANERRORPOLY(1,:));
fprintf(FID, '   %8.4f & %8.4f & %8.4f &  %8.4f & %8.4f  \n',MEANERRORPOLY(2,:));
fclose(FID);

filename=['../results/table4/CompressionPOLYDG_' casenamemaster '_' aggrshock '.tex'];
FID = fopen(filename, 'w');
fprintf(FID, '   %8.0f & %8.0f & %8.0f &  %8.0f & %8.0f & %8.0f \n',Compression(1,:));
% fprintf(FID, '   %8.0f & %8.0f & %8.0f &  %8.0f & %8.0f & %8.0f \n',Compression(2,:));
fclose(FID);

filename=['../results/table4/CompressionPOLYCOEFF_' casenamemaster '_' aggrshock '.tex'];
FID = fopen(filename, 'w');
% fprintf(FID, '   %8.0f & %8.0f & %8.0f &  %8.0f & %8.0f & %8.0f \n',Compression(1,:));
fprintf(FID, '   %8.0f & %8.0f & %8.0f &  %8.0f & %8.0f \n',Compression(2,:));
fclose(FID);
%% Produce results for uncertainty shock
disp('Solve KS model with uncertainty shocks')

casenamemaster='SS_BASELINE_HANC_UNC';

% Select aggregate shock
aggrshock           = 'Uncertainty';
par.rhoS            = 0.84;    % Persistence of variance
par.sigmaS          = 0.54;    % STD of variance shocks

%% Solve for Steady state
disp('Solving Steady State by EGM')
tic
% Set parameters
defineSS_pars_UNC

mainskript_steadystate
toc
%%
disp('Reiter method with State Space Reduction')
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

casename=casenamemaster;

save(['IRF_' casename],'IRF_K','IRF_I','IRF_S','IRF_Y','IRF_C','IRF_Q','IRF_N','IRF_R','IRF_H')

% plot_IRFs

close all
