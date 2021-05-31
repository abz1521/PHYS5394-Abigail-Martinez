x = [2,0.3,10];

a1 = [1,10];
a2 = [-5,0.5];
a3 = [-50,50];

step_size = 0.3;
snr = 10;

samplFreq = 1024; %Hz
nSamples = 2048;
dataLen = nSamples/samplFreq;
% time samples vector
dataX = (0:(nSamples-1))/samplFreq;

% psd function handle
tPSD = @(f) (f>=0 & f<=200).*(200-f).*((sin(f/6)*0.3+1))/10000 + 1;
% set up frequency
kNyq = floor(nSamples/2)+1;
posFreq = (0:(kNyq-1))*(1/dataLen); % vector of positive frequencies
% psd vector
psdVec = tPSD(posFreq);

% generate colored noise
filtOrdr = 100;
rng('default');
noiseVec = statgaussnoisegen(nSamples,[posFreq(:),psdVec(:)],filtOrdr,samplFreq);

% generate qc and normalize to snr
sigVec = crcbgenqcsig(dataX,snr,x);
% norm of the signal squared is inner product of signal with itself
normSigSqrd = innerprodpsd(sigVec,sigVec,samplFreq,psdVec);
sigVec = snr*sigVec/sqrt(normSigSqrd);

% data vector
dataVec = noiseVec + sigVec;


%% Generate matrices
A = a1(1):step_size:a1(2);
nRows = length(A);

% Column of num1 values
num1 = (A-a1(1))./(a1(2)-a1(1));
col1 = num1';
% Columns of same num2 and num3 value respectively
num2 = (x(2)-a2(1))/(a2(2)-a2(1));
col2 = zeros(length(num1),1) + num2;
num3 = (x(3)-a3(1))/(a3(2)-a3(1));
col3 = zeros(length(num1),1) + num3;


mm = [col1,col2,col3];% master matrix

%% Fitness values of mm
P = struct('dataX',dataX,...
                'dataXSq',dataX.^2,...
                'dataXCb',dataX.^3,...
                'dataY',dataVec,...
                'samplFreq',samplFreq,...
                'psdVec',psdVec,...
                'snr',snr,...
                'rmin',[a1(1),a2(1),a3(1)],...
                'rmax',[a1(2),a2(2),a3(2)]);

glrts = glrtqcsig4pso(mm,P);

idmin = find(glrts == min(glrts));

figure;
plot(A,glrts,'-p','MarkerIndices',[idmin],'MarkerFaceColor','red','MarkerSize',10);
title('Fitness Values vs QC Parameter Range');
xlabel('Parameter Value a1');
ylabel('Fitness Value (-GLRT)');
legend({['Minimum at a1 = ',num2str(A(idmin))]},'Location','South East');
