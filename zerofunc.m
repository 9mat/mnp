function [ y, d ] = zerofunc( x )
%ZEROFUNC Summary of this function goes here
%   Detailed explanation goes here

y = 0;

if nargout > 1
    d = zeros(size(x));
end

end

