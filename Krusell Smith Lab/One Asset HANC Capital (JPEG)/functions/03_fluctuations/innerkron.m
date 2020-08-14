function [ R ] = innerkron( n_f,n_v,fmat, varargin )
%innerkron Calculates expressions such as fmat(A*B*C....), where fmat is a
%matrix with symmetric columns. A, B, C, ..., are matrices, and * is a
%Kronecker product. These expressions arise in high order multivariate
%chain rules of the composite function f(v(x)). In this case, fmat is an
%array of high derivatives of f w.r.t v, reshaped into a matrix with the
%same number of rows as f. A, B, C, etc would be derivatives of v w.r.t x.
%
% © Copyright, Oren Levintal, June 13, 2016.

R=fmat;

for j=1:length(varargin)
    R=reshape(R,numel(R)/n_v,n_v);
    R=(R*reshape(varargin{j},n_v,numel(varargin{j})/n_v))';
end

R=reshape(R,numel(R)/n_f,n_f)';

end

