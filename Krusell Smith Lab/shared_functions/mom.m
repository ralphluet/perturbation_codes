function [sigyJ,sigxJ]=mom(gx,hx,varshock,J, method)
%[sigyJ,sigxJ]=mom(gx,hx,varshock,J, method)
% Computes the unconditional variance-covariance matrix of x(t) with x(t+J), that is sigxJ=E[x(t)*x(t+J)'], 
%and the unconditional variance covariaance matrix of y(t) with y(t+J), that is sigyJ=E[y(t)*y(t+J)']
% where x(t) evolves as
% x(t+1) = hx x(t) + e(t+1)
%and y(t) evolves according to 
% y(t) = gx x(t)
%where Ee(t)e(t)'=varshock
%The parameter J can be any integer
%method =1 : use doubling algorithm
%method neq 1 : use algebraic method
%(c) Stephanie Schmitt-Grohe and Martin Uribe, April 18, 1990, renewed January 24, 2000 and August 18, 2001. 


if nargin<4
J=0;
end


if nargin<5
    method =1;
end


if method == 1 
%disp('method=doubling')

%Doubling algorithm
hx_old=hx;
sig_old=varshock;
sigx_old=eye(size(hx));
diferenz=.1;
while diferenz>1e-25;
sigx=hx_old*sigx_old*hx_old'+sig_old;

diferenz = max(max(abs(sigx-sigx_old)));
sig_old=hx_old*sig_old*hx_old'+sig_old;
hx_old=hx_old*hx_old;
sigx_old=sigx;
end    %while diferenz


else


%Algebraic method
%Get the variance of x
disp('method=kronecker')


sigx = zeros(size(hx));
F = kron(hx,hx);
sigx(:) = (eye(size(F))-F)\varshock(:);


end   %if method



%Get E{x(t)*x(t+J)'}
sigxJ=hx^(-min(0,J))*sigx*(hx')^(max(0,J));


%Get E{y(t)*y(t+J)'}
sigyJ=real(gx*sigxJ*gx');
sigxJ=real(sigxJ);