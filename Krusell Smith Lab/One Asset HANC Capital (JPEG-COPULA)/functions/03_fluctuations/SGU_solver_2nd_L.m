function [gxx_gss,hxx_hss,gss,hss,F11,F12,F13,F14,F21,F22,F23,F24,F31,F32,F33,F34,F41,F42,F43,F44,par] = SGU_solver_2nd_L(F,F1,F2,F3,F4,gx,hx,numstates,numcontrols,oc,mpar,par,grid,eta0)

% numstates=mpar.numstates; numcontrols=mpar.numcontrols;
%now do numerical linearization
State       = zeros(numstates,1);
State_m     = State;
Contr       = zeros(numcontrols,1);
Contr_m     = Contr;
tic
[Fb,~,~]  = F(State,State_m,Contr,Contr_m);%Fsys(State,State,Contr,Contr,Xss,Yss,Gamma_state,Gamma_control,InvGamma,par,mpar,grid,targets);
toc
%Use Schmitt Gohe Uribe Algorithm
% E[x' u']' =inv(A)*B*[x u]'
% A = [dF/dx' dF/du'], B =[dF/dx dF/du]
% A = [F1 F2]; B=[F3 F4]

p = gcp('nocreate');
pool='local'; % select machine

if isempty(p)
    c = parcluster;
    parpool(pool,c.NumWorkers-1)
    p = gcp('nocreate');
end

p = gcp('nocreate');


F11 = zeros(numstates+numcontrols,numstates,numstates);
F12 = zeros(numstates+numcontrols,numstates,numcontrols);
F13 = zeros(numstates+numcontrols,numstates,numstates);
F14 = zeros(numstates+numcontrols,numstates,numcontrols);

F21 = zeros(numstates+numcontrols,numcontrols,numstates);
F22 = zeros(numstates+numcontrols,numcontrols,numcontrols);
F23 = zeros(numstates+numcontrols,numcontrols,numstates);
F24 = zeros(numstates+numcontrols,numcontrols,numcontrols);

F31 = zeros(numstates+numcontrols,numstates,numstates);
F32 = zeros(numstates+numcontrols,numstates,numcontrols);
F33 = zeros(numstates+numcontrols,numstates,numstates);
F34 = zeros(numstates+numcontrols,numstates,numcontrols);

F41 = zeros(numstates+numcontrols,numcontrols,numstates);
F42 = zeros(numstates+numcontrols,numcontrols,numcontrols);
F43 = zeros(numstates+numcontrols,numcontrols,numstates);
F44 = zeros(numstates+numcontrols,numcontrols,numcontrols);


% parameters which control numerical differentiation and Reiter solution
numscale = 1; %number of step sizes for numerical differentiation
pnum     = p.NumWorkers;

packagesize = ceil(numstates/(3*pnum));
blocks      = ceil(numstates/packagesize);

% Absolute deviations
par.scaleval1 = 1e-5; %vector of numerical differentiation step sizes
par.scaleval2 = 1e-5; %vector of numerical differentiation step sizes

tic
disp('Computing Hessian F11=D^2F/DXprime^2')
disp(['Total number of parallel blocks: ' num2str(blocks) '.'])

parfor bl=1:blocks
    range = ((bl-1)*packagesize+1):min(packagesize*bl,numstates);
    DF11 = zeros(length(Fb),length(range),numstates);
    cc = zeros(numcontrols,1);
    ss = zeros(numstates,1);
    for Xpct1 = range
        for Xpct2 = 1:numstates
        Xp1 = zeros(numstates,1); mXp1 = zeros(numstates,1);
        Xp2 = zeros(numstates,1); mXp2 = zeros(numstates,1);
        h1 = par.scaleval1;
        Xp1(Xpct1) = h1; mXp1(Xpct1) = -h1;
        Xp2(Xpct2) = h1; mXp2(Xpct2) = -h1;
        F_xp1xp2 = F(Xp1+Xp2,ss,cc,cc);
        F_xp1mxp2 = F(Xp1+mXp2,ss,cc,cc);
        F_mxp1xp2 = F(mXp1+Xp2,ss,cc,cc);
        F_mxp1mxp2 = F(mXp1+mXp2,ss,cc,cc);
        DF11(:,Xpct1-(bl-1)*packagesize,Xpct2) = (F_xp1xp2 - F_xp1mxp2 - F_mxp1xp2 + F_mxp1mxp2)/(4*h1^2);
        end
    end
   
    
    FF11{bl}=DF11;
    disp(['Block number: ' num2str(bl) ' done.'])
end

for i=1:ceil(numstates/packagesize)
    range = ((i-1)*packagesize+1):min(packagesize*i,numstates);
    F11(:,range,:)=FF11{i};
end
toc


tic
disp('Computing Hessian F12=D^2F/DXpDYp')
disp(['Total number of parallel blocks: ' num2str(blocks) '.'])

parfor bl=1:blocks
    range = ((bl-1)*packagesize+1):min(packagesize*bl,numstates);
    DF12 = zeros(length(Fb),length(range),numcontrols);
    cc = zeros(numcontrols,1);
    ss = zeros(numstates,1);
    for Xpct = range
        for Ypct = 1:numcontrols
        Xp = zeros(numstates,1); mXp = zeros(numstates,1);
        h1 = par.scaleval1;
        Xp(Xpct) = h1; mXp(Xpct) = -h1;
        Yp  = zeros(numcontrols,1); mYp = zeros(numcontrols,1);
        h2 = par.scaleval1;
        Yp(Ypct) = h2; mYp(Ypct) = -h2;
        F_xp1yp1 = F(Xp,ss,Yp,cc);
        F_xp1myp1 = F(Xp,ss,mYp,cc);
        F_mxp1yp1 = F(mXp,ss,Yp,cc);
        F_mxp1myp1 = F(mXp,ss,mYp,cc);
        DF12(:,Xpct-(bl-1)*packagesize,Ypct) = (F_xp1yp1 - F_xp1myp1 - F_mxp1yp1 + F_mxp1myp1)/(4*h1*h2);
        end
    end
  
    FF12{bl}=DF12;

    disp(['Block number: ' num2str(bl) ' done.'])
end

for i=1:ceil(numstates/packagesize)
    range = ((i-1)*packagesize+1):min(packagesize*i,numstates);
    F12(:,range,:)=FF12{i};
end
toc


tic
disp('Computing Hessian F13=D^2F/DXprimeDX')
disp(['Total number of parallel blocks: ' num2str(blocks) '.'])

parfor bl=1:blocks
    range = ((bl-1)*packagesize+1):min(packagesize*bl,numstates);
    DF13 = zeros(length(Fb),length(range),numstates);
    cc = zeros(numcontrols,1);
    ss = zeros(numstates,1);
    for Xpct = range
        for Xct = 1:numstates
        Xp = zeros(numstates,1); mXp = zeros(numstates,1);
        h1 = par.scaleval1;
        Xp(Xpct) = h1; mXp(Xpct) = -h1;
        X  = zeros(numstates,1); mX = zeros(numstates,1);
        h2 = par.scaleval1;
        X(Xct) = h2; mX(Xct) = -h2;
        F_xp1x1 = F(Xp,X,cc,cc);
        F_xp1mx1 = F(Xp,mX,cc,cc);
        F_mxp1x1 = F(mXp,X,cc,cc);
        F_mxp1mx1 = F(mXp,mX,cc,cc);
        
        DF13(:,Xpct-(bl-1)*packagesize,Xct) = (F_xp1x1 - F_xp1mx1 - F_mxp1x1 + F_mxp1mx1)/(4*h1*h2);
        end
    end

    FF13{bl}=DF13;

    disp(['Block number: ' num2str(bl) ' done.'])
end

for i=1:ceil(numstates/packagesize)
    range = ((i-1)*packagesize+1):min(packagesize*i,numstates);
    F13(:,range,:)=FF13{i};
end
toc


tic
disp('Computing Hessian F14=D^2F/DXpDY')
disp(['Total number of parallel blocks: ' num2str(blocks) '.'])

parfor bl=1:blocks
    range = ((bl-1)*packagesize+1):min(packagesize*bl,numstates);
    DF14 = zeros(length(Fb),length(range),numcontrols);
    cc = zeros(numcontrols,1);
    ss = zeros(numstates,1);
    for Xpct = range
        for Yct = 1:numcontrols
        Xp = zeros(numstates,1); mXp = zeros(numstates,1);
        h1 = par.scaleval1;
        Xp(Xpct) = h1; mXp(Xpct) = -h1;
        Y  = zeros(numcontrols,1); mY = zeros(numcontrols,1);
        h2 = par.scaleval1;
        Y(Yct) = h2; mY(Yct) = -h2;

        F_xp1y1 = F(Xp,ss,cc,Y);
        F_xp1my1 = F(Xp,ss,cc,mY);
        F_mxp1y1 = F(mXp,ss,cc,Y);
        F_mxp1my1 = F(mXp,ss,cc,mY);

        DF14(:,Xpct-(bl-1)*packagesize,Yct) = (F_xp1y1 - F_xp1my1 - F_mxp1y1 + F_mxp1my1)/(4*h1*h2);
        end
    end
   
    FF14{bl}=DF14;

    disp(['Block number: ' num2str(bl) ' done.'])
end

for i=1:ceil(numstates/packagesize)
    range = ((i-1)*packagesize+1):min(packagesize*i,numstates);
    F14(:,range,:)=FF14{i};
end
toc

packagesize = ceil(numcontrols/(3*pnum));
blocks=ceil(numcontrols/packagesize);
%%
tic
disp('Computing Hessian F21=D^2F/DYpDXp')
disp(['Total number of parallel blocks: ' num2str(blocks) '.'])

parfor bl=1:blocks
    range = ((bl-1)*packagesize+1):min(packagesize*bl,numcontrols);
    DF21 = zeros(length(Fb),length(range),numstates);
    cc = zeros(numcontrols,1);
    ss = zeros(numstates,1);
    for Ypct = range
        for Xpct = 1:numstates
        Yp = zeros(numcontrols,1); mYp = zeros(numcontrols,1);
        h1 = par.scaleval1;
        Yp(Ypct) = h1; mYp(Ypct) = -h1;
        Xp  = zeros(numstates,1); mXp = zeros(numstates,1);
        h2 = par.scaleval1;
        Xp(Xpct) = h2; mXp(Xpct) = -h2;
        F_yp1xp1 = F(Xp,ss,Yp,cc);
        F_yp1mxp1 = F(mXp,ss,Yp,cc);
        F_myp1xp1 = F(Xp,ss,mYp,cc);
        F_myp1mxp1 = F(mXp,ss,mYp,cc);
        DF21(:,Ypct-(bl-1)*packagesize,Xpct) = (F_yp1xp1 - F_yp1mxp1 - F_myp1xp1 + F_myp1mxp1)/(4*h1*h2);
        end
    end
    
    FF21{bl}=DF21;
    disp(['Block number: ' num2str(bl) ' done.'])
end

for i=1:ceil(numcontrols/packagesize)
    range = ((i-1)*packagesize+1):min(packagesize*i,numcontrols);
    F21(:,range,:)=FF21{i};
end
toc

%%
tic
disp('Computing Hessian F22=D^2F/DXp^2')
disp(['Total number of parallel blocks: ' num2str(blocks) '.'])

parfor bl=1:blocks
    range = ((bl-1)*packagesize+1):min(packagesize*bl,numcontrols);
    DF22 = zeros(length(Fb),length(range),numcontrols);
    cc = zeros(numcontrols,1);
    ss = zeros(numstates,1);
    for Ypct1 = range
        for Ypct2 = 1:numcontrols
        Yp1 = zeros(numcontrols,1); mYp1 = zeros(numcontrols,1);
        Yp2 = zeros(numcontrols,1); mYp2 = zeros(numcontrols,1);
        h1 = par.scaleval1;
        Yp1(Ypct1) = h1; mYp1(Ypct1) = -h1;
        Yp2(Ypct2) = h1; mYp2(Ypct2) = -h1;
        F_yp1yp2 = F(ss,ss,Yp1+Yp2,cc);
        F_yp1myp2 = F(ss,ss,Yp1+mYp2,cc);
        F_myp1yp2 = F(ss,ss,mYp1+Yp2,cc);
        F_myp1myp2 = F(ss,ss,mYp1+mYp2,cc);
        DF22(:,Ypct1-(bl-1)*packagesize,Ypct2) = (F_yp1yp2 - F_yp1myp2 - F_myp1yp2 + F_myp1myp2)/(4*h1^2);
        end
    end
    
    FF22{bl}=DF22;

    disp(['Block number: ' num2str(bl) ' done.'])
end

for i=1:ceil(numcontrols/packagesize)
    range = ((i-1)*packagesize+1):min(packagesize*i,numcontrols);
    F22(:,range,:)=FF22{i};
end
toc


%%
tic
disp('Computing Hessian F23=D^2F/DYpDX')
disp(['Total number of parallel blocks: ' num2str(blocks) '.'])

parfor bl=1:blocks
    range = ((bl-1)*packagesize+1):min(packagesize*bl,numcontrols);
    DF23 = zeros(length(Fb),length(range),numstates);
    cc = zeros(numcontrols,1);
    ss = zeros(numstates,1);
    for Ypct = range
        for Xct = 1:numstates
        Yp = zeros(numcontrols,1); mYp = zeros(numcontrols,1);
        h1 = par.scaleval1;
        Yp(Ypct) = h1; mYp(Ypct) = -h1;
        X  = zeros(numstates,1); mX = zeros(numstates,1);
        h2 = par.scaleval1;
        X(Xct) = h2; mX(Xct) = -h2;
        F_yp1x1 = F(ss,X,Yp,cc);
        F_yp1mx1 = F(ss,mX,Yp,cc);
        F_myp1x1 = F(ss,X,mYp,cc);
        F_myp1mx1 = F(ss,mX,mYp,cc);
        DF23(:,Ypct-(bl-1)*packagesize,Xct) = (F_yp1x1 - F_yp1mx1 - F_myp1x1 + F_myp1mx1)/(4*h1*h2);
        end
    end
    
    FF23{bl}=DF23;
    disp(['Block number: ' num2str(bl) ' done.'])
end

for i=1:ceil(numcontrols/packagesize)
    range = ((i-1)*packagesize+1):min(packagesize*i,numcontrols);
    F23(:,range,:)=FF23{i};
end
toc

%%

tic
disp('Computing Hessian F24=D^2F/DYpDY')
disp(['Total number of parallel blocks: ' num2str(blocks) '.'])

parfor bl=1:blocks
    range = ((bl-1)*packagesize+1):min(packagesize*bl,numcontrols);
    DF24 = zeros(length(Fb),length(range),numcontrols);
    cc = zeros(numcontrols,1);
    ss = zeros(numstates,1);
    for Ypct = range
        for Yct = 1:numcontrols
        Yp = zeros(numcontrols,1); mYp = zeros(numcontrols,1);
        h1 = par.scaleval1;
        Yp(Ypct) = h1; mYp(Ypct) = -h1;
        Y  = zeros(numcontrols,1); mY = zeros(numcontrols,1);
        h2 = par.scaleval1;
        Y(Yct) = h2; mY(Yct) = -h2;
        F_yp1y1 = F(ss,ss,Yp,Y);
        F_yp1my1 = F(ss,ss,Yp,mY);
        F_myp1y1 = F(ss,ss,mYp,Y);
        F_myp1my1 = F(ss,ss,mYp,mY);
        DF24(:,Ypct-(bl-1)*packagesize,Yct) = (F_yp1y1 - F_yp1my1 - F_myp1y1 + F_myp1my1)/(4*h1*h2);
        end
    end
 
    FF24{bl}=DF24;

    disp(['Block number: ' num2str(bl) ' done.'])
end

for i=1:ceil(numcontrols/packagesize)
    range = ((i-1)*packagesize+1):min(packagesize*i,numcontrols);
    F24(:,range,:)=FF24{i};
end
toc


packagesize = ceil(numstates/(3*pnum));
blocks      = ceil(numstates/packagesize);

tic
disp('Computing Hessian F31=D^2F/DXDXp')
disp(['Total number of parallel blocks: ' num2str(blocks) '.'])

parfor bl=1:blocks
    range = ((bl-1)*packagesize+1):min(packagesize*bl,numstates);
    DF31 = zeros(length(Fb),length(range),numstates);
    cc = zeros(numcontrols,1);
    ss = zeros(numstates,1);
    for Xct = range
        for Xpct = 1:numstates
        X = zeros(numstates,1); mX = zeros(numstates,1);
        h1 = par.scaleval1;
        X(Xct) = h1; mX(Xct) = -h1;
        Xp  = zeros(numstates,1); mXp = zeros(numstates,1);
        h2 = par.scaleval1;
        Xp(Xpct) = h2; mXp(Xpct) = -h2;
        F_x1xp1 = F(Xp,X,cc,cc);
        F_x1mxp1 = F(mXp,X,cc,cc);
        F_mx1xp1 = F(Xp,mX,cc,cc);
        F_mx1mxp1 = F(mXp,mX,cc,cc);
        
        DF31(:,Xct-(bl-1)*packagesize,Xpct) = (F_x1xp1 - F_x1mxp1 - F_mx1xp1 + F_mx1mxp1)/(4*h1*h2);
        end
    end
    FF31{bl}=DF31;

    disp(['Block number: ' num2str(bl) ' done.'])
end

for i=1:ceil(numstates/packagesize)
    range = ((i-1)*packagesize+1):min(packagesize*i,numstates);
    F31(:,range,:)=FF31{i};
end
toc

%%
tic
disp('Computing Hessian F32=D^2F/DXDYp')
disp(['Total number of parallel blocks: ' num2str(blocks) '.'])

parfor bl=1:blocks
    range = ((bl-1)*packagesize+1):min(packagesize*bl,numstates);
    DF32 = zeros(length(Fb),length(range),numcontrols);
    cc = zeros(numcontrols,1);
    ss = zeros(numstates,1);
    for Xct = range
        for Ypct = 1:numcontrols
        X = zeros(numstates,1); mX = zeros(numstates,1);
        h1 = par.scaleval1;
        X(Xct) = h1; mX(Xct) = -h1;
        Yp  = zeros(numcontrols,1); mYp = zeros(numcontrols,1);
        h2 = par.scaleval1;
        Yp(Ypct) = h2; mYp(Ypct) = -h2;
        F_x1yp1 = F(ss,X,Yp,cc);
        F_x1myp1 = F(ss,X,mYp,cc);
        F_mx1yp1 = F(ss,mX,Yp,cc);
        F_mx1myp1 = F(ss,mX,mYp,cc);
        DF32(:,Xct-(bl-1)*packagesize,Ypct) = (F_x1yp1 - F_x1myp1 - F_mx1yp1 + F_mx1myp1)/(4*h1*h2);
        end
    end
    
    FF32{bl}=DF32;

    disp(['Block number: ' num2str(bl) ' done.'])
end

for i=1:ceil(numstates/packagesize)
    range = ((i-1)*packagesize+1):min(packagesize*i,numstates);
    F32(:,range,:)=FF32{i};
end
toc

%%
tic
disp('Computing Hessian F33=D^2F/DX^2')
disp(['Total number of parallel blocks: ' num2str(blocks) '.'])

parfor bl=1:blocks
    range = ((bl-1)*packagesize+1):min(packagesize*bl,numstates);
    DF33 = zeros(length(Fb),length(range),numstates);
    cc = zeros(numcontrols,1);
    ss = zeros(numstates,1);
    for Xct1 = range
        for Xct2 = 1:numstates
        X1 = zeros(numstates,1); mX1 = zeros(numstates,1);
        X2 = zeros(numstates,1); mX2 = zeros(numstates,1);
        h1 = par.scaleval1;
        X1(Xct1) = h1; mX1(Xct1) = -h1;
        X2(Xct2) = h1; mX2(Xct2) = -h1;
        F_x1x2 = F(ss,X1+X2,cc,cc);
        F_x1mx2 = F(ss,X1+mX2,cc,cc);
        F_mx1x2 = F(ss,mX1+X2,cc,cc);
        F_mx1mx2 = F(ss,mX1+mX2,cc,cc);
        DF33(:,Xct1-(bl-1)*packagesize,Xct2) = (F_x1x2 - F_x1mx2 - F_mx1x2 + F_mx1mx2)/(4*h1^2);
        end
    end
 
    FF33{bl}=DF33;
    disp(['Block number: ' num2str(bl) ' done.'])
end

for i=1:ceil(numstates/packagesize)
    range = ((i-1)*packagesize+1):min(packagesize*i,numstates);
    F33(:,range,:)=FF33{i};

end
toc



%%
tic
disp('Computing Hessian F34=D^2F/DXDY')
disp(['Total number of parallel blocks: ' num2str(blocks) '.'])

parfor bl=1:blocks
    range = ((bl-1)*packagesize+1):min(packagesize*bl,numstates);
    DF34 = zeros(length(Fb),length(range),numcontrols);
    cc = zeros(numcontrols,1);
    ss = zeros(numstates,1);
    for Xct = range
        for Yct = 1:numcontrols
        X = zeros(numstates,1); mX = zeros(numstates,1);
        h1 = par.scaleval1;
        X(Xct) = h1; mX(Xct) = -h1;
        Y  = zeros(numcontrols,1); mY = zeros(numcontrols,1);
        h2 = par.scaleval1;
        Y(Yct) = h2; mY(Yct) = -h2;

        F_x1y1 = F(ss,X,cc,Y);
        F_x1my1 = F(ss,X,cc,mY);
        F_mx1y1 = F(ss,mX,cc,Y);
        F_mx1my1 = F(ss,mX,cc,mY);

        DF34(:,Xct-(bl-1)*packagesize,Yct) = (F_x1y1 - F_x1my1 - F_mx1y1 + F_mx1my1)/(4*h1*h2);
        end
    end
    
    FF34{bl}=DF34;

    disp(['Block number: ' num2str(bl) ' done.'])
end

for i=1:ceil(numstates/packagesize)
    range = ((i-1)*packagesize+1):min(packagesize*i,numstates);
    F34(:,range,:)=FF34{i};
end
toc


packagesize = ceil(numcontrols/(3*pnum));
blocks=ceil(numcontrols/packagesize);
%%
tic
disp('Computing Hessian F41=D^2F/DYDXp')
disp(['Total number of parallel blocks: ' num2str(blocks) '.'])

parfor bl=1:blocks
    range = ((bl-1)*packagesize+1):min(packagesize*bl,numcontrols);
    DF41 = zeros(length(Fb),length(range),numstates);
    cc = zeros(numcontrols,1);
    ss = zeros(numstates,1);
    for Yct = range
        for Xpct = 1:numstates
        Y = zeros(numcontrols,1); mY = zeros(numcontrols,1);
        h1 = par.scaleval1;
        Y(Yct) = h1; mY(Yct) = -h1;
        Xp  = zeros(numstates,1); mXp = zeros(numstates,1);
        h2 = par.scaleval1;
        Xp(Xpct) = h2; mXp(Xpct) = -h2;
        F_y1xp1 = F(Xp,ss,cc,Y);
        F_y1mxp1 = F(mXp,ss,cc,Y);
        F_my1xp1 = F(Xp,ss,cc,mY);
        F_my1mxp1 = F(mXp,ss,cc,mY);
        DF41(:,Yct-(bl-1)*packagesize,Xpct) = (F_y1xp1 - F_y1mxp1 - F_my1xp1 + F_my1mxp1)/(4*h1*h2);
        end
    end
    FF41{bl}=DF41;
    disp(['Block number: ' num2str(bl) ' done.'])
end

for i=1:ceil(numcontrols/packagesize)
    range = ((i-1)*packagesize+1):min(packagesize*i,numcontrols);
    F41(:,range,:)=FF41{i};
end
toc

%%

tic
disp('Computing Hessian F42=D^2F/DYDYp')
disp(['Total number of parallel blocks: ' num2str(blocks) '.'])

parfor bl=1:blocks
    range = ((bl-1)*packagesize+1):min(packagesize*bl,numcontrols);
    DF42 = zeros(length(Fb),length(range),numcontrols);
    cc = zeros(numcontrols,1);
    ss = zeros(numstates,1);
    for Yct = range
        for Ypct = 1:numcontrols
        Y = zeros(numcontrols,1); mY = zeros(numcontrols,1);
        h1 = par.scaleval1;
        Y(Yct) = h1; mY(Yct) = -h1;
        Yp  = zeros(numcontrols,1); mYp = zeros(numcontrols,1);
        h2 = par.scaleval1;
        Yp(Ypct) = h2; mYp(Ypct) = -h2;
        F_y1yp1 = F(ss,ss,Yp,Y);
        F_y1myp1 = F(ss,ss,mYp,Y);
        F_my1yp1 = F(ss,ss,Yp,mY);
        F_my1myp1 = F(ss,ss,mYp,mY);
        DF42(:,Yct-(bl-1)*packagesize,Ypct) = (F_y1yp1 - F_y1myp1 - F_my1yp1 + F_my1myp1)/(4*h1*h2);
        end
    end
 
    FF42{bl}=DF42;

    disp(['Block number: ' num2str(bl) ' done.'])
end

for i=1:ceil(numcontrols/packagesize)
    range = ((i-1)*packagesize+1):min(packagesize*i,numcontrols);
    F42(:,range,:)=FF42{i};
end
toc

%%
tic
disp('Computing Hessian F43=D^2F/DYDX')
disp(['Total number of parallel blocks: ' num2str(blocks) '.'])

parfor bl=1:blocks
    range = ((bl-1)*packagesize+1):min(packagesize*bl,numcontrols);
    DF43 = zeros(length(Fb),length(range),numstates);
    cc = zeros(numcontrols,1);
    ss = zeros(numstates,1);
    for Yct = range
        for Xct = 1:numstates
        Y = zeros(numcontrols,1); mY = zeros(numcontrols,1);
        h1 = par.scaleval1;
        Y(Yct) = h1; mY(Yct) = -h1;
        X  = zeros(numstates,1); mX = zeros(numstates,1);
        h2 = par.scaleval1;
        X(Xct) = h2; mX(Xct) = -h2;
        F_y1x1 = F(ss,X,cc,Y);
        F_y1mx1 = F(ss,mX,cc,Y);
        F_my1x1 = F(ss,X,cc,mY);
        F_my1mx1 = F(ss,mX,cc,mY);
        DF43(:,Yct-(bl-1)*packagesize,Xct) = (F_y1x1 - F_y1mx1 - F_my1x1 + F_my1mx1)/(4*h1*h2);
        end
    end
    
    FF43{bl}=DF43;
    disp(['Block number: ' num2str(bl) ' done.'])
end

for i=1:ceil(numcontrols/packagesize)
    range = ((i-1)*packagesize+1):min(packagesize*i,numcontrols);
    F43(:,range,:)=FF43{i};
end
toc

%%
tic
disp('Computing Hessian F44=D^2F/DY^2')
disp(['Total number of parallel blocks: ' num2str(blocks) '.'])

parfor bl=1:blocks
    range = ((bl-1)*packagesize+1):min(packagesize*bl,numcontrols);
    DF44 = zeros(length(Fb),length(range),numcontrols);
    cc = zeros(numcontrols,1);
    ss = zeros(numstates,1);
    for Yct1 = range
        for Yct2 = 1:numcontrols
        Y1 = zeros(numcontrols,1); mY1 = zeros(numcontrols,1);
        Y2 = zeros(numcontrols,1); mY2 = zeros(numcontrols,1);
        h1 = par.scaleval1;
        Y1(Yct1) = h1; mY1(Yct1) = -h1;
        Y2(Yct2) = h1; mY2(Yct2) = -h1;
        F_y1y2 = F(ss,ss,cc,Y1+Y2);
        F_y1my2 = F(ss,ss,cc,Y1+mY2);
        F_my1y2 = F(ss,ss,cc,mY1+Y2);
        F_my1my2 = F(ss,ss,cc,mY1+mY2);
        DF44(:,Yct1-(bl-1)*packagesize,Yct2) = (F_y1y2 - F_y1my2 - F_my1y2 + F_my1my2)/(4*h1^2);
        end
    end
    
    FF44{bl}=DF44;

    disp(['Block number: ' num2str(bl) ' done.'])
end

for i=1:ceil(numcontrols/packagesize)
    range = ((i-1)*packagesize+1):min(packagesize*i,numcontrols);
    F44(:,range,:)=FF44{i};
end
toc


%% gxx_hxx and gss_hss

[gxx_gss,hxx_hss]=gxx_hxx_gss_hss(F3,F1,F4,F2,F22,F24,F21,F23,F42,F44,F41,F43,F12,F14,F11,F13,F32,F34,F31,F33,hx,gx,eta0);

gss = gxx_gss(:,end);
hss = hxx_hss(:,end);

end

