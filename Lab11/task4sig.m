%% Calculate significance by Computing GLRT of M data realizations under H0

% number of data realizations to compute
m = 50000;
A = 10;
a1=10;
a2=3;
a3=3;

data1 = load('data1.txt','-ascii').'; 
data2 = load('data2.txt','-ascii').';
data3 = load('data3.txt','-ascii').';

nSamples = length(data1);
% Define sampling frequency 
samplFreq = 1024;
dataLen = nSamples/samplFreq;

%% Generate quadratic chirp signal vector
% Time vector
x = (0:(nSamples-1))/samplFreq;
parVec = [10,3,3];

sigVec = crcbgenqcsig(x,A,[a1,a2,a3]);

noisePSD = @(f) (f>=100 & f<=300).*(f-100).*(300-f)/10000 + 1;

% Vector of positive frequencies
kNyq = floor(nSamples/2)+1;
posFreq = (0:(kNyq-1))*(1/dataLen);
psdPosFreq = noisePSD(posFreq);

%% GLRT of M realizations
% glrtH0 vector of the GLRT values for m noise realizations
glrtH0 = zeros(1,m);

% Manually create fir2 filter once to save time
sqrtPSD = sqrt(psdPosFreq);
b = fir2(100,posFreq/(samplFreq/2),sqrtPSD);

for r = 1:m % m realizations
    inNoise = randn(1,nSamples);
    noiseVec = sqrt(samplFreq)*fftfilt(b,inNoise);
    glrtH0(r) = glrtqcsig(noiseVec,samplFreq,psdPosFreq,parVec);
end

% Observed GLRT
obsGamma1 = glrtqcsig(data1, samplFreq, psdPosFreq, parVec);
obsGamma2 = glrtqcsig(data2, samplFreq, psdPosFreq, parVec);
obsGamma3 = glrtqcsig(data3, samplFreq, psdPosFreq, parVec);

% Significance
alpha1 = sum(glrtH0>=obsGamma1)/m;
alpha2 = sum(glrtH0>=obsGamma2)/m;
alpha3 = sum(glrtH0>=obsGamma3)/m;

% Display results
disp([num2str(m),' noise realizations compared']);
disp(['data1 significance = ',num2str(alpha1)]);
disp(['data2 significance = ',num2str(alpha2)]);
disp(['data3 significance = ',num2str(alpha3)]);