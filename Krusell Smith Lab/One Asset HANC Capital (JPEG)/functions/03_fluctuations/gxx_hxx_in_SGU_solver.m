
function [gxx,hxx] = gxx_hxx_in_SGU_solver(F3,F1,F4,F2,F22,F24,F21,F23,F42,F44,F41,F43,F12,F14,F11,F13,F32,F34,F31,F33,hx,gx)

m=0;
nx = size(hx,1); %rows of hx and hxx
ny = size(gx,1); %rows of gx and gxx
n = nx + ny; %length of f
ngxx = nx^2*ny; %elements of gxx

sg = [ny nx nx]; %size of gxx
sh = [nx nx nx]; %size of hxx

Q = zeros(n*nx*(nx+1)/2,n*nx*nx);
q = zeros(n*nx*(nx+1)/2,1);
gxx=zeros(sg);
hxx=zeros(sh);
GXX=zeros(sg);
HXX=zeros(sh);


for i=1:n
for j=1:nx
%for k=1:nx
for k=1:j
m = m+1;

%First Term
q(m,1) = ( shiftdim(F22(i,:,:),1) * gx * hx(:,k) + shiftdim(F24(i,:,:),1) * gx(:,k) + shiftdim(F21(i,:,:),1) *  hx(:,k) + shiftdim(F23(i,:,k),1) )' * gx * hx(:,j); 

% Second term

GXX(:) = kron(ones(nx^2,1),F2(i,:)');

pGXX = permute(GXX,[2 3 1]);
pGXX(:) = pGXX(:) .* kron(ones(nx*ny,1),hx(:,j));
GXX=ipermute(pGXX,[2 3 1]);

pGXX = permute(GXX,[3 1 2]);
pGXX(:) = pGXX(:) .* kron(ones(nx*ny,1),hx(:,k));
GXX=ipermute(pGXX,[3 1 2]);

Q(m,1:ngxx)=GXX(:)';

GXX=0*GXX;

%Third term

HXX(:,j,k) = (F2(i,:) * gx)';

Q(m,ngxx+1:end)=HXX(:)';

HXX = 0*HXX;

%Fourth Term
q(m,1) = q(m,1) + ( shiftdim(F42(i,:,:),1) * gx * hx(:,k) +  shiftdim(F44(i,:,:),1) * gx(:,k) + shiftdim(F41(i,:,:),1) * hx(:,k) +  shiftdim(F43(i,:,k),1) )' * gx(:,j); 

%Fifth Term

GXX(:,j,k)=F4(i,:)';

Q(m,1:ngxx) = Q(m,1:ngxx) + GXX(:)';

GXX = 0*GXX;

%Sixth term
q(m,1) = q(m,1) + ( shiftdim(F12(i,:,:),1) * gx * hx(:,k) + shiftdim(F14(i,:,:),1) * gx(:,k) + shiftdim(F11(i,:,:),1) * hx(:,k) + F13(i,:,k)')' * hx(:,j);


%Seventh Term

HXX(:,j,k)=F1(i,:)';

Q(m,ngxx+1:end) = Q(m,ngxx+1:end) + HXX(:)';

HXX = 0*HXX;

%Eighth Term
q(m,1) = q(m,1) +  shiftdim(F32(i,j,:),1) * gx * hx(:,k) +  shiftdim(F34(i,j,:),1) * gx(:,k) +  shiftdim(F31(i,j,:),1) * hx(:,k) + F33(i,j,k);

end %k 
end %j
end %i

A = temp(nx,ny);

Qt = Q* A;

xt = -Qt\q;
x = A* xt;

gxx(:)=x(1:ngxx);
hxx(:) = x(ngxx+1:end);
end


function A = temp(nx,ny) %subfunction
%function A = temp(nx,ny)
%This function creates a matrix A, size n*nx^2 by n*nx*(nx+1)/2, such that x = A xtilde, 
%where x is the vector  [gxx(:); hxx(:)] and xtilde is a vector containing
%all the elements  gxx(i,j,k) and hxx(i,j,k), respectively, such that j<=
%k, that is, it is a subset of the elements in the vector x that appear only once.
%This is so because gxx(i,j,k) is symmetric with respect to j and k, that is, gxx(i,j,k)=gxx(i,k,j).
%The reason we use this program is that the matrix we invert (which in the earlier versioni of this program used to be
%Q and now is Qt) is of  size n*nx*(nx+1)/2 by n*nx*(nx+1)/2 rather than of size
%n*nx^2 by n*nx^2. This reduces computation time by about 30 percent. (Feb 16, 2004)

Ahxx=zeros(nx^3, nx^2*(nx+1)/2);
Agxx=zeros(ny*nx^2, ny*nx*(nx+1)/2);
mx=0;
my=0;
for k=1:nx;
    for j=k:nx;
        for i=1:nx; 
            mx=mx+1;
            Ahxx((j-1)*nx+i+(k-1)*nx*nx, mx) = 1;
            Ahxx((k-1)*nx+i+(j-1)*nx*nx, mx) = 1;
        end %i=1:nx
        for i=1:ny; 
            my=my+1;
            Agxx((j-1)*ny+i+(k-1)*ny*nx, my) = 1;
            Agxx((k-1)*ny+i+(j-1)*ny*nx, my) = 1;
        end %i=1:ny
    end %j
end %k

A = zeros((nx+ny)*nx^2,(nx+ny)*nx*(nx+1)/2);
A(1:ny*nx^2, 1:ny*nx*(nx+1)/2)=Agxx;
A(ny*nx^2+1:end, ny*nx*(nx+1)/2+1:end)=Ahxx;
end