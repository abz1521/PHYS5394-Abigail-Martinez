% This function generates a normally distributed pseudo-random number
% R = customrandn(M,S), where R is the normally distributed pseudo-random
% number, m is the mean of the normal distribution, and s is the standard
% deviation of the normal distribution.

% Custom normal PDF N(x;mean,std) of Y = std*X+mean, where X has a PDF N(x;0,1).
function randValueN = customrandn(m,s)
randValueN = s*randn() + m;