function [freqHz] = MatsuokaEstimation(obj,tau,T,b,a)
% this function calculate the frequency based on the Matsuoka estimation
% from the 1985 paper

if ( ((tau+T)*b-(tau*a)) / (tau*a) ) < 0
    freqHz = NaN; % ignore non exciting periods
else
    freqHz = (1/(T*2*pi))*sqrt(((tau+T)*b-(tau*a)) / (tau*a)); % in [Hz]
end

end

