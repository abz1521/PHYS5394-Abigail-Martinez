% Generate a function that uniformly distributes a pseudo-random number
% R = customran(a,b). This generates a pseudo-random number R between the
% interval of [a, b], in which a and b are the boundaries of the uniform probability 
% distribution function. 

%Custom uniform PDF U(x;a,b) of Y = (b-a)X+a, where X has a PDF U(x;0,1)

function randValue = customrand(a,b)
randValue = (b-a)*rand() + a;

