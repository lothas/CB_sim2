function [errMSE,R_squar] = NN_perf(obj,targets,outputs,dispFlag,GraphicsFlag,trainOrTest )
% calculate the R^2 and MSE of a fitting neural net

% trainOrTest- either 'train' or 'test'

% calc R^2
err = targets-outputs;
errVar = var(err,0,2);
inputVar = var(targets,0,2);

% R^2 = 1 - (error variance/input variance)
R_squar = 1-(errVar/inputVar);
errMSE = immse(outputs,targets);

if dispFlag
    disp('NN performance:');
    disp(['the R^2 is: ',num2str(R_squar)]);
    disp(['the MSE is: ',num2str(errMSE)]);
end

if GraphicsFlag
    figure;
    plotregression(targets,outputs,trainOrTest);
end

end