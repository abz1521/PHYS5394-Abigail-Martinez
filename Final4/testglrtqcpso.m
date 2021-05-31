%% Generate data realization and test glrtqcpso.m
% define testing parameters
x = [10,3,3]; % qc coefficient vector x
snr = 10;
nRuns = 8;
nSteps = 1000;
% set the [min,max] values for parameter ranges
a1 = [1,180];
a2 = [1,10];
a3 = [1,10];


%% Generate data ralization
% Set up time
samplFreq = 512;%hz
nSamples = 512;
dataLen = nSamples/samplFreq;
% time samples vector
dataX = (0:(nSamples-1))/samplFreq;


noisePSD = @(f) (f>=50 & f<=100).*(f-50).*(100-f)/625 + 1;
% set up frequency
kNyq = floor(nSamples/2)+1;
posFreq = (0:(kNyq-1))*(1/dataLen); % vector of positive frequencies
% psd vector
psdVec = noisePSD(posFreq);
% plot(posFreq,psdVec);

% generate colored noise
filtOrdr = 100;
flatness = floor(nSamples/10+1); % usual length of anomaly for filt ordr 100
noiseVec = statgaussnoisegen(nSamples+flatness,[posFreq(:),psdVec(:)],filtOrdr,samplFreq);
noiseVec = noiseVec(flatness+1:end);

% generate qc and normalize to snr
sigVec = crcbgenqcsig(dataX,1,x);
% norm of the signal squared is inner product of signal with itself
normSigSqrd = innerprodpsd(sigVec,sigVec,samplFreq,psdVec);
sigVec = snr*sigVec/sqrt(normSigSqrd);

dataVec = noiseVec + sigVec;

%% Call glrtqcpso
inParams = struct('dataX',dataX,...
                'dataXSq',dataX.^2,...
                'dataXCb',dataX.^3,...
                'dataY',dataVec,...
                'samplFreq',samplFreq,...
                'psdVec',psdVec,...
                'snr',snr,...
                'rmin',[a1(1),a2(1),a3(1)],...
                'rmax',[a1(2),a2(2),a3(2)]);

psoParams = struct('steps',nSteps);
outStruct = glrtqcpso(inParams,psoParams,nRuns);

%% Presentation

figure;
hold on;
scatter(dataX,dataVec,'Marker','.');
plot(dataX,sigVec,'Color','Red','LineWidth',2.0);
plot(dataX,outStruct.bestSig,'LineStyle','-','LineWidth',4.0, 'Color',[76,153,0]/255)
% figure;
% hold on;
for lpruns = 1:nRuns
    plot(dataX,outStruct.allRunsOutput(lpruns).estSig,'Color',[51,floor(255*lpruns/8),153]/255,'LineWidth',1.0);
end
legend('Data','Signal','Best Estimated Signal',...
                      ['Estimated Signals from ',num2str(nRuns),' runs'],...
                      '','','','','','','');
title('Data, Signal, and Best Estimates');
ylabel('Amplitude');
xlabel('Time (s)');
disp(['Estimated parameters: a1=',num2str(outStruct.bestQcCoefs(1)),...
                          '; a2=',num2str(outStruct.bestQcCoefs(2)),...
                          '; a3=',num2str(outStruct.bestQcCoefs(3))]);