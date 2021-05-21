% Plot Exponentially Damped Sinusoid Signal

% Parameters

snr = 10;
ta = 0.3;
f = 20;
tau = 0.2;
phi = 0;
l = 0.6;

%Sample Intervals

sampleFreq = 10*f;
sampleInterval=1/sampleFreq;
x=0:sampleInterval:1.0;

%Run Function

sigVec = genEDSSig(x,snr,ta,f,tau,phi,l);

%Plot Function

plot(x,sigVec);
title("Exponentially Damped Sinusoid Signal");