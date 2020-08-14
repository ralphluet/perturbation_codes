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

invutil = @(u)(((1-par.xi).*u).^(1/(1-par.xi)));
invmutil = @(mu)((1./mu).^(1/par.xi));


Xss=[squeeze(sum(joint_distr,2)); ... % marginal distribution liquid
    squeeze(sum(joint_distr,1)'); ... % marginal distribution productivity
    0];

Yss=[c_guess(:); log(par.Q); log(targets.Y);...
    log(par.W); log(par.R); log(grid.K)];

os = length(Xss) - (mpar.nm+mpar.nh);
oc = length(Yss) - (mpar.nm*mpar.nh);
%%
Gamma_state=zeros(mpar.nm+mpar.nh,mpar.nm+mpar.nh-2);
for j=1:mpar.nm-1
    Gamma_state(1:mpar.nm,j)=-Xss(1:mpar.nm);
    Gamma_state(j,j)=1-Xss(j);
    Gamma_state(j,j)=Gamma_state(j,j) -sum(Gamma_state(1:mpar.nm,j));
end
bb=mpar.nm;
for j=1:mpar.nh-1
    Gamma_state(bb+(1:mpar.nh),bb+j-1)=-Xss(bb+(1:mpar.nh));
    Gamma_state(bb+j,bb-1+j)=1-Xss(bb+j);
    Gamma_state(bb+j,bb-1+j)=Gamma_state(bb+j,bb-1+j) -sum(Gamma_state(bb+(1:mpar.nh),bb-1+j));
end


%%
aux = reshape(c_guess,[mpar.nm, mpar.nh]);
[aux,DC1] = mydct(aux,1); % do dct-transformation
[aux, DC2] = mydct(aux,2); % do dct-transformation
DC=aux(:);
[~,ind] = sort(abs(DC(:)),'descend');

% [1,1 1,2
%  2,1 2,2
%  3,1
%  4,1]
NNKEEP=[0:mpar.nm-1]';
NNKEEP=[NNKEEP NNKEEP+1];

NNIND=[1:mpar.nm]';
NNIND=[NNIND NNIND+1];

IND=sub2ind([mpar.nm mpar.nh],[1:mpar.nm 1:mpar.nm],[ones(1,mpar.nm) ones(1,mpar.nm)+1]);
IND=reshape(IND,[mpar.nm mpar.nh]);

KEEP=NNKEEP<=DEGREEN;

indexMUdct=IND(KEEP)';

% indexMUdct=1:keepnumber;

%%
aux= size(Gamma_state); %used for distributions

mpar.numstates   = aux(2)+os;
mpar.numcontrols = length(indexMUdct)+oc;
keepnumber = length(indexMUdct);
State       = zeros(mpar.numstates,1);
State_m     = State;
Contr       = zeros(mpar.numcontrols,1);
Contr_m     = Contr;
%%

F = @(a,b,c,d)Fsys(a,b,c,d,Xss,Yss,Gamma_state,indexMUdct,DC1,DC2,Copula,par,mpar,grid,meshes,P_H,aggrshock,oc,os);

[Fss,LHS,RHS,Distr] = F(State,State_m,Contr,Contr_m);


%% Solve RE via Schmitt-Grohe-Uribe Form
% disp('Take Numerical Derivatives and Solve RE via Schmitt-Grohe-Uribe Form')

[hx,gx,F1,F2,F3,F4,par] = SGU_solver(F, mpar.numstates,mpar.numcontrols,oc,mpar,par,p);

