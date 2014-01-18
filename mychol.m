function [ L ] = mychol( A )
%MYCHOL Summary of this function goes here
%   Detailed explanation goes here

if size(A,1) == 1
    L = sqrt(A(1,1));
    return;
end

[L, q] = chol(A, 'lower');

if q == 0
    return
end

if q == 1
    L = mychol(A(2:end, 2:end));
    L = padarray(L, [1,1], 0 , 'pre');
    return;
end

R = L\A(1:q-1,q:end);
B = mychol(A(q:end, q:end) - R'*R);
L = [L, zeros(size(R)); R', B];

end

