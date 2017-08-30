function [MSE,R2,slope] = NN_perf_calc(obj,targ,outputs,plotFlag)
% calculates the MSE, R^2, Slope

% calc MSE:
MSE = immse(targ,outputs);

% calc R^2
err = targ - outputs;
errVar = var(err,0,2);
inputVar = var(targ,0,2);
% % % R^2 = 1 - (error variance/input variance) % % %
R2 = 1-(errVar/inputVar);

% calc slope:
coefficients = polyfit(targ, outputs, 1);
slope = coefficients(1);

%% plotting stuff:
if plotFlag
    disp('NN performance:');
    disp(['the R^2 is: ',num2str(R2)]);
    disp(['the MSE is: ',num2str(MSE)]);
    disp(['the slope is: ',num2str(slope)]);
    
    figure;
    plotregression(targ,outputs,'out Vs targ');
end

end


