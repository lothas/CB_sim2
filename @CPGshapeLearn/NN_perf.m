function [errMSE,R_squar] = NN_perf(obj,varargin)
% calculate the R^2 and MSE of a fitting neural net

% trainOrTest- either 'train' or 'test'
switch nargin
    case 6
        targets = varargin{1};
        outputs = varargin{2};
        dispFlag = varargin{3};
        GraphicsFlag = varargin{4};
        trainOrTest = varargin{5};
    case 4
        targets = varargin{1};
        outputs = varargin{2};
        dispFlag = varargin{3};
        GraphicsFlag = false;
        trainOrTest = 'test';
    case 3
        targets = varargin{1};
        outputs = varargin{2};
        dispFlag = false;
        GraphicsFlag = false;
    otherwise
        disp('invalid number of inputs...')
        error(['the inputs order is: 1)targets, 2)outputs',...
        '3)dispFlag, 4)GraphicsFlag','trainOrTest']);
end

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