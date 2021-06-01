function [fitVal,varargout] = glrtqcsig4pso(xVec,params)
[nVecs,~]=size(xVec);

%storage for fitness values
fitVal = zeros(nVecs,1);

%Check for out of bound coordinates and flag them
validPts = crcbchkstdsrchrng(xVec);
%Set fitness for invalid points to infty
fitVal(~validPts)=inf;
xVec(validPts,:) = s2rv(xVec(validPts,:),params);

for lpc = 1:nVecs
    if validPts(lpc)
    % Only the body of this block should be replaced for different fitness
    % functions
        x = xVec(lpc,:);
        fitVal(lpc) = ssrqc(x, params);
    end
end

%Return real coordinates if requested
if nargout > 1
    varargout{1}=xVec;
end

%LLR after maximizing over amplitude parameter
function ssrVal = ssrqc(x,params)
%Generate normalized quadratic chirp
phaseVec = x(1)*params.dataX + x(2)*params.dataXSq + x(3)*params.dataXCb;
qc = sin(2*pi*phaseVec);
% Inner product of signal with itself
normSigSqrd = innerprodpsd(qc,qc,params.samplFreq,params.psdVec);
% Normalization factor
normFac = params.snr/sqrt(normSigSqrd);
% Normalize qc signal
qc = normFac*qc;

%Compute fitness
ssrVal = -innerprodpsd(params.dataY,qc,params.samplFreq,params.psdVec)^2;