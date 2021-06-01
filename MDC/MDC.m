% read TRAINIGDATA.MAT and ANALYSISDATA.MAT
trainingmat = load('TrainingData.mat');
dataT = trainingmat.trainData;
dataVec = load('analysisData.mat').dataVec;

% Parameter Range
a1 = [40,100];
a2 = [1,50];
a3 = [1,15];

% gen time vector
nSamples = 2048;
sampFreq = 1024; % Hz
dataLen = nSamples/sampFreq;
dataX = (0:(nSamples-1))/sampFreq;

% estimate PSD
[pxx,posFreq] = pwelch(dataT, 1024,[],[],sampFreq);
%***************
% interpolate psd and freq vec to get length 1024
psdVec = interp1(1:length(pxx),pxx,linspace(1,length(pxx),1025),'cubic');
posFreq = interp1(1:length(posFreq),posFreq,linspace(1,length(posFreq),1025),'linear');
sqrtPSD = sqrt(psdVec);
b = fir2(100,posFreq/(sampFreq/2),sqrtPSD);
%% Calculate the GLRT
% m noise realizations
m = 1000;
glrtH0 = zeros(1,m);
pars = [70,25,7];

for r = 1:m % m realizations
    inNoise = randn(1,nSamples);
    noiseVec = sqrt(sampFreq)*fftfilt(b,inNoise);
    glrtH0(r) = glrtqcsig(noiseVec,sampFreq,psdVec,pars(:));
end

[mSign,mIndex] = min(glrtH0);
disp(['Lowest significance=',num2str(mSign),newline,...
'Parameters a1=',num2str(pars(1)),...
         '; a2=',num2str(pars(2)),...
         '; a3=',num2str(pars(3))]);
% after testing several different parameters,
% I believe there's a signal bc of very low significance     

%% Call PSO
% parameter structrue
nRuns = 8;
nSteps = 2000;
snr = 1;
inParams = struct('dataX',dataX,...
                  'dataXSq',dataX.^2,...
                  'dataXCb',dataX.^3,...
                  'dataY',dataVec,...
                  'samplFreq',sampFreq,...
                  'psdVec',psdVec,...
                  'snr',snr,...
                  'rmin',[a1(1),a2(1),a3(1)],...
                  'rmax',[a1(2),a2(2),a3(2)]);

outStruct = glrtqcpso(inParams,struct('maxSteps',2000),nRuns);

%% Estimation of SNR
normSigSqrd = innerprodpsd(outStruct.bestSig,outStruct.bestSig,sampFreq,psdVec);
% Normalize best signal to specified SNR
sigVec = snr*outStruct.bestSig/sqrt(normSigSqrd);

% Obtain LLR values under H0 for bestSig
flatness = 100;
nH0Data = 500;
llrH0 = zeros(1,nH0Data);
for lp = 1:nH0Data
    inNoise = randn(1,nSamples+flatness);
    noiseVec = sqrt(sampFreq)*fftfilt(b,inNoise);
    noiseVec = noiseVec(flatness+1:end);
    llrH0(lp) = innerprodpsd(noiseVec,sigVec,sampFreq,psdVec);
end

llrH1 = innerprodpsd(dataVec,sigVec,sampFreq,psdVec);

%Signal to noise ratio estimate
estSNR = (llrH1-mean(llrH0))/std(llrH0);
bestSig = estSNR*sigVec;

%% Plot
figure;
hold on;
plot(dataX,dataVec,'Color','Magenta');
% plot(dataX,sigVec,'Color','Red','LineWidth',2.0);
plot(dataX,bestSig,'LineStyle','-','LineWidth',1.0, 'Color',[76,153,0]/255);
title('Data and Estimated Signal');
xlabel('Time (s)');
ylabel('Amplitude');
legend('Data','Estimated Signal');
disp(['PSO estimated parameters: a1=',num2str(outStruct.bestQcCoefs(1)),...
                              '; a2=',num2str(outStruct.bestQcCoefs(2)),...
                              '; a3=',num2str(outStruct.bestQcCoefs(3)),...
                              newline,'Estimated SNR=',num2str(estSNR)]);