function sigVec = genEDSSig(x,snr,ta,f,tau,phi,l)
% Generate an Exponentially Damped Sinusoid Signal
% X is the vector of the time stamps at which the value of the signal is to
% be computed

%SNR * e^-(x-ta)sin(2*pi*f*t)

% Abigail Martinez, Jan 25, 2021

sigVec=snr*exp(-(x-ta)/tau).*sin(2*pi*f*x + phi);

for i=1:length(x)
    if x(i) < ta
        sigVec(i)=0;
    elseif x(i) > ta+l
        sigVec(i)=0;
    end
end
