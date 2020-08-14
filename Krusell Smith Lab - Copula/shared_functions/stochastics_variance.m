function [P_H,grid,par]=stochastics_variance(par, mpar, grid)
% stochastics_variance generates transition probabilities

% Values for transition and productivities from JEDC comparison project
P_H=[0.6 0.4; 0.044445 0.955555];

grid.h=[par.nu 1];

Paux=P_H^1000^1000;
hh=Paux(1,1:mpar.nh)*grid.h(1:mpar.nh)';
par.H=hh(1);
par.N = par.H;

end