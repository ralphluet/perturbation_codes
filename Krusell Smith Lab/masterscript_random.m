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
    %     c.NumWorkers
    parpool(pool,10)
    p = gcp('nocreate');
end

%% Select aggregate shock
aggrshock           = 'TFP';
par.rhoS            = 0.75;     % Persistence
par.sigmaS          = 0.007;    % STD
mpar.ns             = 3;        % number of TFP states
mpar.T(1,1)         = 1500;     % number of periods
mpar.T_KS=mpar.T(1,1);
mpar.T(2,1)         = 500;      % burn-in phase
mpar.maxlag         = 40;       % number of periods IRF

mpar.nK             = 3;        % number of grid points for KS


%% Parameters

% Household Parameters
par.beta        = 0.99; %0.95 - 0.995     % Discount factor
par.xi          = 1;   %1-4       % CRRA
par.gamma       = 1;   % 0.5 - 2       % Inverse Frisch elasticity

% Income Process
par.rhoH        = 0.9; % 0.7 - 0.95   % Persistence of productivity
par.sigmaH      = 0.2; % 0.05 - 0.4   % STD of productivity shocks

% Firm Side Parameters
par.alpha       = 0.64;% 0.5 - 0.75       % Labor share 2/3
par.delta       = 0.1/4;     % Depreciation rate
par.phi         = 0;        % Capital adj costs

% Price of Capital in SS
par.Q  = 1;

% Labor supply
par.nu= 0.15; %productivity when "UE", i.e. income when UE


% Grids
% Idiosyncratic States
mpar.nm         = 100;
mpar.nk         = 0;
mpar.nh         = 2;
mpar.tauchen    ='importance';


% Numerical Parameters
mpar.crit    = 1e-10;


NoP=3;

% Household Parameters
% par.beta        = 0.99; %0.95 - 0.995     % Discount factor
beta=linspace(0.95,0.99,NoP);
% par.xi          = 1;   %1-4       % CRRA
xi=linspace(1,4,NoP);
% par.gamma       = 1;   % 0.5 - 2       % Inverse Frisch elasticity
gamma=linspace(0.5,2,NoP);

% Income Process
% par.rhoH        = 0.9; % 0.7 - 0.95   % Persistence of productivity
rhoH=linspace(0.7,0.95,NoP);
% par.sigmaH      = 0.2; % 0.05 - 0.4   % STD of productivity shocks
sigmaH=linspace(0.05,0.4,NoP);

% Firm Side Parameters
% par.alpha       = 0.64;% 0.5 - 0.75       % Labor share 2/3
alpha=linspace(0.5,0.75,NoP);

% Shocks
rhoS=linspace(0.5,0.95,NoP);
sigmaS=linspace(0.001,0.02,NoP);

%%

count=0;


for NN1=1:NoP
    for NN2=1:NoP
        for NN3=1:NoP
            for NN4=1:NoP
                for NN5=1:NoP
                    for NN6=1:NoP
                        for NN7=1:NoP
                            for NN8=1:NoP
                                
                                count=count+1;
                                
                                disp(count)
                                disp(count/NoP^8)
                                
                                par.beta=beta(NN1);
                                par.xi=xi(NN2);
                                par.gamma=gamma(NN3);
                                par.rhoH=rhoH(NN4);
                                par.sigmaH=sigmaH(NN5);
                                par.alpha=alpha(NN6);
                                par.rhoS = rhoS(NN7);
                                par.sigmaS = sigmaS(NN8);
                                
       
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
                                
                                
                                DIST_MEAN_RR(count,:)=mean(abs([RD.^0.5]));
                                DIST_MAX_RR(count,:)=max(abs([RD.^0.5]));
                                
                                %%
                                disp('Krusell Smith method')
                                cd('One Asset HANC Capital (KS)')
                                
                                tic
                                mainskript
                                toc
                                cd ..
                                addpath(genpath('One Asset HANC Capital (KS)/functions'))
                                
                                [IRF_K_KS, IRF_S_KS]=KSsimulationIRF(k_star,P_H,joint_distr,grid,mesh,mpar,par);
                                
                                [Kmat(:,3)]=KSsimulation(k_star,P_H,PS_S,pr_s_master(:),joint_distr,grid,mesh,mpar);
                                Kaux=Kmat(:,3);
                                Imat(:,3)=(Kaux(2:end)-(1-par.delta)*Kaux(1:end-1));
                                
                                Klinear_KS=zeros(mpar.T(1,1),1);
                                Klinear_KS(1)=targets.K;
                                for t=1:mpar.T(1,1)-1
                                    
                                    Klinear_KS(t+1)=exp(par.beta_K(1,smaster(t))   + par.beta_K(2,smaster(t)).*log(Klinear_KS(t)));
                                    
                                    
                                end
                                
                                DenHaanError(:,3)=100*log(Kaux./Klinear_KS);
                                
                                DH_MEAN_RR(count,:)= mean(abs(DenHaanError));
                                DH_MAX_RR(count,:)= max(abs(DenHaanError));
                                
                                path(oldpath)
                                                                
                            end
                        end
                    end
                end
            end
        end
    end
end

save PARADH_ALL count DH_MEAN_RR DH_MAX_RR
%%
filename=['../results/table10/METAAvgAbsErrordenHaan_' casenamemaster '_' aggrshock '.tex'];
FID = fopen(filename, 'w');
fprintf(FID, '   %8.4f & %8.4f & %8.4f \n', mean(DH_MEAN_RR));
fclose(FID);
filename=['../results/table10/METAMaxAbsErrordenHaan_' casenamemaster '_' aggrshock '.tex'];
FID = fopen(filename, 'w');
fprintf(FID, '   %8.4f & %8.4f & %8.4f \n',mean(DH_MAX_RR));
fclose(FID);





%%
CC=reshape(DH_MAX_RR,[3,3,3^6,3]);
mCC=mean(CC,3);
maxCC=max(CC,[],3)-mCC;
minCC=min(CC,[],3)-mCC;
figurename=['MAX_DH_Comparison_' casenamemaster  '_' aggrshock];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 500 1000])
set(gca, 'FontName','arial','FontSize',40);
set(gca,'Position',[0.125 0.15 0.825 0.825])
for j = 1:3
    subplot(3,1,j), grouped_error_bar(squeeze(mCC(j,:,:,:)),squeeze(maxCC(j,:,:,:)),squeeze(minCC(j,:,:,:)))
    title(['\sigma_S =' num2str(sigmaS(j))]) 
    if j==1; legend({'Reiter-Reduction','Reiter-Full', 'Krusell-Smith','Min-Max range difference across experiments'}, 'Location','NorthWest'); end
    set(gca,'XTickLabel',{'\rho_S = 0.5','\rho_S = 0.725','\rho_S = 0.95'});
end
printpdf(gcf,['../results/figure8/' figurename])
cccc=maxCC-minCC;



CC=reshape(DH_MEAN_RR,[3,3,3^6,3]);
mCC=mean(CC,3);
maxCC=max(CC,[],3)-mCC;
minCC=min(CC,[],3)-mCC;
figurename=['MEAN_DH_Comparison_' casenamemaster  '_' aggrshock];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 500 1000])
set(gca, 'FontName','arial','FontSize',40);
set(gca,'Position',[0.125 0.15 0.825 0.825])
for j = 1:3
    subplot(3,1,j), grouped_error_bar(squeeze(mCC(j,:,:,:)),squeeze(maxCC(j,:,:,:)),squeeze(minCC(j,:,:,:)))
    title(['\sigma_S =' num2str(sigmaS(j))])
    if j==1; legend({'Reiter-Reduction','Reiter-Full', 'Krusell-Smith','Min-Max range across experiments'}, 'Location','NorthWest'); end
    set(gca,'XTickLabel',{'\rho_S = 0.5','\rho_S = 0.725','\rho_S = 0.95'});
end
printpdf(gcf,['../results/figure8/' figurename])
cccc=maxCC-minCC;
function grouped_error_bar(M,maxS,minS)
    ngroups = size(M, 1);
    nbars = size(M, 2);
    % Calculating the width for each bar group
    groupwidth = min(0.8, nbars/(nbars + 1.5));
    hold on
    bar(M,'grouped')
    xticks(1:ngroups)
    for i = 1:nbars
        x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
        errorbar(x, M(:,i), minS(:,i),maxS(:,i), 'k.','LineWidth',2);
    end
   
    hold off
end

