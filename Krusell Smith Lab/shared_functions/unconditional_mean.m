function [Ey,Ex] = unconditional_mean(gx, hx, gxx, hxx, gss, hss, eta, sig)
%
%[Ey,Ex] = unconditional_mean(gx, hx, gxx, hxx, gss, hss, eta, sig) 
%computes the unconditional mean of the state vector
% x and the co-state vector y
%
%Inputs: Arrays defining the second-order expansion of x and y
% as described in Schmitt-Grohe and Uribe (JEDC, 2004). Namely,
% hx is nx by nx
% gx is ny by nx
% hxx is nx by nx by nx
% gxx is ny by nx by nx
% hss is nx by 1
% gss is ny by 1
% eta is nx by ne
% sig is a parameter scaling the innovation to the exogenous variables 
%
%Output:  vectors Ex (nx by 1) and Ey (ny by 1) containing unconditonal means for x and y 
%
%Calls MOM.M
% (c) Stephanie Schmitt-Grohe and Martin Uribe, April 23, 2004


%Number of shocks
ne = size(eta,2);
%Columns of X and Y
nx=size(hx,1);
ny=size(gx,1);

%Unconditional var/cov matrix of state vector x
[var_cov_y,var_cov_x]=mom(gx,hx,sig^2*(eta*eta'),0); % eta=eta0;

%The unconditional means of x and y can be written as the solutions for Ex and Ey of
%Ex=hx*Ex+a/2
%Ey=gx*Ey+b/2;
%where a and b are given by
for i=1:nx
a(i,1) =  sum(sum( var_cov_x .* squeeze(hxx(i,:,:)) )) + hss(i) *sig^2; 
end %i

for i=1:ny
b(i,1) = sum(sum( var_cov_x .* squeeze(gxx(i,:,:)) )) + gss(i) * sig^2; 
end %i

Ex=(eye(nx)-hx)\(a/2);
Ey=gx*Ex+b/2;
end