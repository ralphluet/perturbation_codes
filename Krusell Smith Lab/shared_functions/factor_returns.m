function [ W_fc,WW,RB,R_fc,Y ] = factor_returns(K,par,mpar)
%factor_returns
  
%%  GHH Preferences    
  
    W_fc =  par.alpha.* (K./par.N).^(1-par.alpha);        
    R_fc = (1-par.alpha).* (par.N./K).^(par.alpha)- par.delta;
    Y = par.N.^(par.alpha).*K.^(1-par.alpha);
 
    %%
    NW=par.gamma/(1+par.gamma).*(par.N/par.H).*W_fc;
    WW=NW*ones(mpar.nm,mpar.nh); %Wages
    RB = (1+R_fc);

end

