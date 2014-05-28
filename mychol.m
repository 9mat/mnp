function [ L ] = mychol( A )
%MYCHOL Modified Choleski decomposition
%   Perform Choleski decomposition even when the matrix is singular or near
%   singular (but still positive semidefinite). 

% Base case: a scalar
if isscalar(A)
    assert(A > -1e-6);
    if A < 0; A = 0; end;
    L = sqrt(A);
    return;
end

% Call Matlab decomposition
[L, q] = chol(A, 'lower');

% Ideal case: strictly positive definite
if q == 0
    return
end

% If the decomposition fails right at first element, the first row and the 
% first column must be closed to zeros; do the decomposition recuresily for
% the rest of the matrix
if q == 1
    assert(mean(abs(A(1,:))) + mean(abs(A(:,1))) < 1e-5);
    L = mychol(A(2:end, 2:end));
    l = size(L,1);
    L = [zeros(l+1,1), [zeros(1,l); L] ];
    return;
end

% Normal case: do the decomposition recursively
R = L\A(1:q-1,q:end);
B = mychol(A(q:end, q:end) - R'*R);
L = [L, zeros(size(R)); R', B];

end

