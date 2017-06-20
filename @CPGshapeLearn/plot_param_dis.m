function plot_param_dis(obj,whichOne)
% plot the parameter disturbution

binNum = 100;
inputsNames = obj.inputsNames;
outputNames = obj.outputNames;

figure;

switch whichOne
    case {'in','Inputs','inputs'}
        graphsNum = length(inputsNames);
        param = inputsNames;
    case {'out','Outputs','outputs'}
        graphsNum = length(outputNames);
        param = outputNames;
    case {'all','everything'}
        param = {'tau','b','a','c','A0','A1','A2','A3','B1','B2','B3','freq'};
        graphsNum = length(param);
    otherwise
        graphsNum = 1;
        param = cell(1,1);
        param{1,1} = whichOne;
end

switch graphsNum
    case 1
        temp = obj.get_param(param{1,1});
        hist(temp,binNum);    grid minor;
        title(['historiogram of ',param{1,1}]);
    otherwise
        for j=1:graphsNum
            subplot(3,ceil(graphsNum/3),j);
            temp = obj.get_param(param{1,j});
            hist(temp,binNum);    grid minor;
            title(['historiogram of ',param{1,j}]);
        end
end

end

