function [fitVal,varargout] = glrtqcsig4pso(xVec,P)

[nVecs,~]=size(xVec);

fitVal = zeros(nVecs,1);

boundary_chk = crcbchkstdsrchrng(xVec);
%Set fitness for invalid points to infty
fitVal(~boundary_chk)=inf;
xVec(boundary_chk,:) = s2rv(xVec(boundary_chk,:),P);

for lpc = 1:nVecs
    if boundary_chk(lpc)
        x = xVec(lpc,:);
        phaseVec = x(1)*P.dataX + x(2)*P.dataXSq + x(3)*P.dataXCb;
        qc = sin(2*pi*phaseVec);
        normSigSqrd = innerprodpsd(qc,qc,P.samplFreq,P.psdVec);
        normFac = P.snr/sqrt(normSigSqrd); % Factor to normalize qc
        qc = normFac*qc; % Normalize qc

        % fitness calc
        fitVal = -innerprodpsd(P.dataY,qc,P.samplFreq,P.psdVec)^2;
    end
end