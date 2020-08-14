function [Y,X] = simu_2nd(gx, hx, gxx, hxx, gss, hss, eta, sig, x0, e)
%
%[Y,X] = simu_2nd(gx, hx, gxx, hxx, gss, hss, eta, sig, x0, e) simulates times series from the model:  
%
%xp^i = hx^i_a x_a + 1/2 hxx^i_ab x_a x_b + 1/2 hss^i sig^2 + sig * eta^i_c ep_c
%
%y^i = gx^i_a x_a + 1/2 gxx^i_ab x_a x_b + 1/2 gss^i sig^2
%
%where 
% hx is nx by nx
% gx is ny by nx
% hxx is nx by nx by nx
% gxx is ny by nx by nx
% hss is nx by 1
% gss is ny by 1
% eta is nx by ne
% sig is a parameter scaling the innovation to the exogenous variables 
% e is T by ne random shock
% T is the length of the simulation 
%
% x0 is the initial condition for X
%
%Output:  vectors X (T by nx) and Y (T by ny) containing time series for x and y and e (T by 1) containing  
%
% (c) Stephanie Schmitt-Grohe and Martin Uribe, January 22, 2002
% eta=eta0; e=e_simul; x0=x0_simul;

%Number of periods
T=size(e,1)+1;
%Number of shocks
ne = size(eta,2);
%Columns of X and Y
nx=size(hx,1);
ny=size(gx,1);

x0=x0(:);

X(1,1:nx)=x0';

%Decompose x into its first- and second-order components
xf=x0;
xs=0*x0;

for t=1:T

%Construct y(t)
   for i=1:ny
     y(i,1) = gx(i,:) * (xf+xs) + 1/2 * xf' * squeeze(gxx(i,:,:)) * xf + 1/2 * gss(i,1)*sig^2; 
   end %i

Y(t,1:ny) = y';

%Remember xf(t) and xs(t) when entering t+1
xoldf = xf; 
xolds = xs;

%Construct x(t+1)
   if t<T
     for i=1:nx
       xf(i,1) = hx(i,:) * xoldf + sig * eta(i,:) * e(t,:)';
       xs(i,1) = hx(i,:) * xolds + 1/2 * xoldf' * squeeze(hxx(i,:,:)) * xoldf + 1/2 * hss(i,1)*sig^2;
     end
   X(t+1,:)=(xf + xs)';
   end %if t<T

end

end