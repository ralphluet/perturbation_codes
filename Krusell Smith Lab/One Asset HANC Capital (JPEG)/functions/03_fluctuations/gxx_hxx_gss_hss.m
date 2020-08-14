function [gxx_gss,hxx_hss]=gxx_hxx_gss_hss(F3,F1,F4,F2,F22,F24,F21,F23,F42,F44,F41,F43,F12,F14,F11,F13,F32,F34,F31,F33,hx,gx,x0)

n_f = size(F4,1);
n_y = size(F4,2);
n_x = size(F3,2)+1;
n_x2 = 0;
n_v = (n_y+n_x)*2;
n_e = 1;

%% fv
efxp = [F1,zeros(n_f,1)];
efx  = [F3,zeros(n_f,1)];
fv = [F2,F4,efxp,efx];


%%
efypxp = cat(3,F21,zeros(n_f,n_y,1));
efypx = cat(3,F23,zeros(n_f,n_y,1));

mfypyp=reshape(F22,[n_f*n_y,n_y]);
mfypy=reshape(F24,[n_f*n_y,n_y]);
mfypxp=reshape(efypxp,[n_f*n_y,n_x]);
mfypx=reshape(efypx,[n_f*n_y,n_x]);
%%
efyxp = cat(3,F41,zeros(n_f,n_y,1));
efyx = cat(3,F43,zeros(n_f,n_y,1));

mfyyp=reshape(F42,[n_f*n_y,n_y]);
mfyy=reshape(F44,[n_f*n_y,n_y]);
mfyxp=reshape(efyxp,[n_f*n_y,n_x]);
mfyx=reshape(efyx,[n_f*n_y,n_x]);

%%
efxpyp = cat(2,F12,zeros(n_f,1,n_y));
efxpy = cat(2,F14,zeros(n_f,1,n_y));
efxpxp = cat(3,cat(2,F11,zeros(n_f,1,n_x-1)),zeros(size(F11,1),n_x,1));
efxpx = cat(3,cat(2,F13,zeros(n_f,1,n_x-1)),zeros(size(F13,1),size(F13,2)+1,1));

mfxpyp=reshape(efxpyp,[n_f*n_x,n_y]);
mfxpy=reshape(efxpy,[n_f*n_x,n_y]);
mfxpxp=reshape(efxpxp,[n_f*n_x,n_x]);
mfxpx=reshape(efxpx,[n_f*n_x,n_x]);

%%
efxyp = cat(2,F32,zeros(size(F32,1),1,size(F32,3)));
efxy = cat(2,F34,zeros(size(F34,1),1,size(F34,3)));
efxxp = cat(3,cat(2,F31,zeros(size(F31,1),1,size(F31,2))),zeros(size(F31,1),size(F31,2)+1,1));
efxx = cat(3,cat(2,F33,zeros(size(F33,1),1,size(F33,2))),zeros(size(F33,1),size(F33,2)+1,1));

mfxyp=reshape(efxyp,[n_f*n_x,n_y]);
mfxy=reshape(efxy,[n_f*n_x,n_y]);
mfxxp=reshape(efxxp,[n_f*n_x,n_x]);
mfxx=reshape(efxx,[n_f*n_x,n_x]);

fvv = [mfypyp,mfypy,mfypxp,mfypx;...
       mfyyp,mfyy,mfyxp,mfyx;...
       mfxpyp,mfxpy,mfxpxp,mfxpx;...
       mfxyp,mfxy,mfxxp,mfxx];

fvv = fvv(:);

gx = [gx,zeros(n_y,1)];
hx = [[hx;zeros(1,n_x-1)],zeros(n_x,1)]; hx(end,end)=1;

M2=1;
eta_ = [x0;0];

eta2_M2=reshape([eta_*reshape(M2,n_e,n_e)]',n_e,n_x);
eta2_M2=reshape([eta_*eta2_M2]',n_x^2,1);

% unique=nchoosek(n_x+2-1,2);
% unique=unique-nchoosek(n_x-1+1-1,1);

Vx0=[gx*hx;gx;hx;speye(n_x)];
Vx1=[gx;sparse(n_y,n_x);speye(n_x,n_x);sparse(n_x,n_x)];

Ezeta2=[ sparse(n_x^2,n_x^2-1) , eta2_M2 ];

A=innerkron(n_f,n_v,fvv,Vx0,Vx0)+innerkron(n_f,n_v,fvv,Vx1,Vx1)*Ezeta2;

fy_fxp_fypgx=[fv(:,n_y+1:2*n_y) fv(:,2*n_y+1:2*n_y+n_x)+fv(:,1:n_y)*gx];

G=fy_fxp_fypgx(:,1:n_f);
H=fy_fxp_fypgx(:,n_f+1:n_f+n_x2);

D=sparse(n_f,n_f);
D(:,1:n_y)=fv(:,1:n_y);

C=A;

%Block 1: xx
spx=sparse([ones(n_x-1,1);0]);
choosex2=kron(spx,spx);
choosex2=logical(choosex2==1);
tempeye=speye(n_x^2);
Z=tempeye(:,choosex2);
CZ=C*Z;
hx_hat=hx(1:end-1,1:end-1);

G=full(G);
D=full(D);
hx_hat=full(hx_hat);

[~,Xtemp]=gensylv( 2,G,D,hx_hat,full(-CZ) );


% acc=norm(full(CZ+AkronkC(D*Xtemp,hx_hat,2)+G*Xtemp));
% if acc>1e-8
% warning(['Second order derivatives may be inaccurate, norm(error)=' num2str(acc) '. Try a different solver.'])
% end

X=zeros(n_f,n_x^2);
X(:,choosex2)=full(Xtemp); 
gxx_hat=Xtemp(1:n_y,:);
hat_eta2_M2=eta2_M2(choosex2,:);

%Block 2: ss
sps=sparse([zeros(n_x-1,1);1]);
choosex2=kron(sps,sps);
choosex2=logical(choosex2==1);
Z=tempeye(:,choosex2);
CZ=C*Z;

Xtemp=-(D+G)\(CZ+(F2*gxx_hat)*hat_eta2_M2);

acc=norm(full(CZ+(F2*gxx_hat)*hat_eta2_M2+(D+G)*Xtemp));
if acc>1e-8
warning(['Second order derivatives may be inaccurate, norm(error)=' num2str(acc) '. Try a different solver.'])
end
Xtemp=reshape(Xtemp,[],size(Z,2));
X(:,choosex2)=full(Xtemp); clear Xtemp

gxx_gss=X(1:n_y,:);
hxx_hss=[X(n_y+1:end,:);zeros(1,n_x^2)];
end
