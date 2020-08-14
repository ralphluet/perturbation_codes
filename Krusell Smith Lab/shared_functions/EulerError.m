









function [EulerE] = EulerError(kprime,Y,R,mutil,par,mpar,gri,prob,meshes)
    % EulerError calculates the error from the Euler Equation for the stocahstic
    % growth model given a savings policy KPRIME (dimensions: k x z).
    % Y (dimensions: k x z) is income for all combinations of assets (k)
    % and productivity (z) on grid GRI.k and GRI.z respectively.
    % PROB is the Transition probability matrix (dimensions: z x z')
    % and MESHES the meshes of gri.z and gri.k.

kprime=reshape(kprime,[mpar.nk,mpar.nz]); % make sure that kprime comes as a matrix
if min(Y(:)-kprime(:))>0 && min(kprime(:))>0 % feasible policy
    MU  = mutil(Y-kprime); %Calculate marginal utility given consumption (y(k,z)-k'(k,z))
    % Calculate expected marginal utility E (mu(k',z')) (dimensions k' x z)
    EMU = par.beta*(R.*MU)*prob.z';
    % use this to define an interpolant
    EMUfun  = griddedInterpolant({gri.k,gri.z},EMU,'linear');
    EMU_val = EMUfun(kprime,meshes.z); % Evaluate interpolant at k'(k,z) to obtain f(k,z) = E_z mu(k'(k,z),z')

    EE = MU-EMU_val; % Write up Euler Equation Error
    EulerE=EE(:); % Stack as vector
else % infeasible policy
    EulerE=ones(size(Y(:)))*1e15;
end
end