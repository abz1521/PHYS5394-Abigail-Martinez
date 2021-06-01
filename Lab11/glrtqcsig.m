function llr = glrtqcsig(dataVec, sampFreq, psdVec, parVec)
%% Parameters for data realization
% Number of samples and sampling frequency.
nSamples = 2048;
sampFreq = 1024;
timeVec = (0:(nSamples-1))/sampFreq;

%% Generate  data realization
% Noise + SNR=10 signal. 
a1=10;
a2=3;
a3=3;
A=10;
sigVec = crcbgenqcsig(timeVec,A,[a1,a2,a3]);
% Signal normalized to SNR=10
% [sig4data,~]=normsig4psd(sig4data,sampFreq,psdPosFreq,10);
% dataVec = noiseVec+sig4data;


%% Compute GLRT
%Generate the unit norm signal (i.e., template). Here, the value used for
%'A' does not matter because we are going to normalize the signal anyway.
%Note: the GLRT here is for the unknown amplitude case, that is all other
%signal parameters are known
sigVec = crcbgenqcsig(timeVec,A,[a1,a2,a3]);
%We do not need the normalization factor, just the  template vector
[templateVec,~] = normsig4psd(sigVec,sampFreq,psdVec,10);
% Calculate inner product of data with template
llr = innerprodpsd(dataVec,templateVec,sampFreq,psdVec);
%GLRT is its square
llr = llr^2;

end

