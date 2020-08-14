
invutil = @(u)(((1-par.xi).*u).^(1/(1-par.xi)));
invmutil = @(mu)((1./mu).^(1/par.xi));


Xss=[squeeze(sum(joint_distr,2)); ... % marginal distribution liquid
    squeeze(sum(joint_distr,1)'); ... % marginal distribution productivity
    (joint_distr(:));...
    0];

Yss=[c_guess(:);  log(par.Q); log(targets.Y);...
    log(par.W); log(par.R); log(grid.K)];

os = length(Xss) - (mpar.nm+mpar.nh+(mpar.nm)*(mpar.nh));
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


%% DCT Transformation

[indexMUdct,DC1,DC2,needed,XXMU] = do_dct(c_guess,mpar.nm,mpar.nh, 0.9999);


[indexCOPdct,DC3,DC4,needed,XXCOP] = do_dct((joint_distr(1:end-1,1:end-1)),mpar.nm-1,mpar.nh-1, 0.95);


%%
aux= size(Gamma_state); %used for distributions

mpar.numstates   = aux(2)+os+length(indexCOPdct);
mpar.numcontrols = length(indexMUdct)+oc;
State       = zeros(mpar.numstates,1);
State_m     = State;
Contr       = zeros(mpar.numcontrols,1);
Contr_m     = Contr;
