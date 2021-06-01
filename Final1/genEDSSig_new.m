function sigVec = genEDSSig_new(dataX,snr,P)

% sigVec=snr*exp(-(dataX-ta)/tau).*sin(2*pi*f*dataX + phi);

phaseVec = 2*pi*P.length*P.freq*dataX+P.phi;
EDSSig = snr*exp(-(dataX-P.ta)/P.tau).*sin(2*pi*P.freq*dataX + P.phi);
EDSSig(dataX<P.ta) = 0;
EDSSig(dataX>P.ta+P.length) = 0;
sinSig = sin(phaseVec);
sigVec = EDSSig.*sinSig;
