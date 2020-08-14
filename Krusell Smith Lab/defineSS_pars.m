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


%% Grids
% Idiosyncratic States
mpar.nm         = 100;
mpar.nk         = 0;
mpar.nh         = 2;
mpar.tauchen    ='importance';


%% Numerical Parameters
mpar.crit    = 1e-10;
