function outResults = glrtqcpso(inParams,psoParams,nRuns)


nSamples = length(inParams.dataX);

fHandle = @(x) glrtqcsig4pso(x,inParams);

nDim = 3;
outStruct = struct('bestLocation',[],...
                   'bestFitness', [],...
                   'totalFuncEvals',[]);
                    
outResults = struct('allRunsOutput',struct('fitVal', [],...
                                           'qcCoefs',zeros(1,3),...
                                           'estSig',zeros(1,nSamples),...
                                           'totalFuncEvals',[]),...
                    'bestRun',[],...
                    'bestFitness',[],...
                    'bestSig', zeros(1,nSamples),...
                    'bestQcCoefs',zeros(1,3));

%Allocate storage for outputs: results from all runs are stored
for lpruns = 1:nRuns
    outStruct(lpruns) = outStruct(1);
    outResults.allRunsOutput(lpruns) = outResults.allRunsOutput(1);
end
%Independent runs of PSO in parallel. Change 'parfor' to 'for' if the
%parallel computing toolbox is not available.
parfor lpruns = 1:nRuns
    %Reset random number generator for each worker
    rng(lpruns);
    outStruct(lpruns)=crcbpso(fHandle,nDim,psoParams);
end

%Prepare output
fitVal = zeros(1,nRuns);
for lpruns = 1:nRuns   
    fitVal(lpruns) = outStruct(lpruns).bestFitness;
    outResults.allRunsOutput(lpruns).fitVal = fitVal(lpruns);
    [~,qcCoefs] = fHandle(outStruct(lpruns).bestLocation);
    outResults.allRunsOutput(lpruns).qcCoefs = qcCoefs;
    Sig = crcbgenqcsig(inParams.dataX,1,qcCoefs);
    % Normalize estimated qc signal (colored noise)
    normSigSq = innerprodpsd(Sig,inParams.samplFreq,inParams.psdVec);
    % Normalization factor (estimated amplitude)
    normFac = 1/sqrt(normSigSq);
    Sig = normFac*Sig;
    estAmp = innerprodpsd(inParams.dataY,Sig,inParams.samplFreq,inParams.psdVec);
    Sig = estAmp*Sig;
    outResults.allRunsOutput(lpruns).estSig = Sig;
    outResults.allRunsOutput(lpruns).totalFuncEvals = outStruct(lpruns).totalFuncEvals;
end
%Find the best run
[~,bestRun] = min(fitVal(:));
outResults.bestRun = bestRun;
outResults.bestFitness = outResults.allRunsOutput(bestRun).fitVal;
outResults.bestSig = outResults.allRunsOutput(bestRun).estSig;
outResults.bestQcCoefs = outResults.allRunsOutput(bestRun).qcCoefs;