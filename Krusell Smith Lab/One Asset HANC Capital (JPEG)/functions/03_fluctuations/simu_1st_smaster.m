function [Y,X, e] = simu_1st_smaster(gx, hx, eta, T,x0, e, smaster)
%
% [Y,X, e] = simu_1st(gx, hx, eta, T,x0, e);
%simulates times series from the model:  
%
%x_t+1 = hx x_t + eta e_t+1 for t=1,...,T-1
%
%y_t = gx * x_t for t=1,...,T
%
%where 
% hx is nx by nx
% gx is ny by nx
% eta is nx by ne
% e is T by ne random shock (default~N(0,I)). Note that e(1) is not used in the simulation because X(1) is taken as given.
% T is the length of the simulation, default 500
%
% x0 is the initial condition for X, default 0
%
%Output:  vectors X, (T by nx), Y (T by ny) 
%containing time series for x and y and e (T 
%by ne) containing the innovations  
%
% (c) Stephanie Schmitt-Grohe and Martin Uribe, December 12, 2007

nx=size(hx,1); %number of states
ny=size(gx,1); %number of controls
ne = size(eta,2); %Number of shocks

%Set default initial state
if nargin <5
x0 = zeros(nx,1);
end

%Setdefault sample size
if nargin<4
T = 500;
end

if nargin<6
e = randn(T,ne);
end


%initialize X and Y
X = zeros(T, nx);
Y = zeros(T, ny); 

x0=x0(:);



X(1,1:nx)=x0';
Y(1,1:ny) = X(1,:)*gx';

for t=2:T
X(t,:) = X(t-1,:)*hx';
Y(t,1:ny) = X(t,:)*gx';
end

end