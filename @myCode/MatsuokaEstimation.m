function [freqHz] = MatsuokaEstimation(tau,T,b,a)
% this function calculate the frequency based on the Matsuoka estimation
% from the 1985 paper

freqHz = (1/(T*2*pi))*sqrt(((tau+T)*b-(tau*a))/(tau*a)); % in [Hz]

end

