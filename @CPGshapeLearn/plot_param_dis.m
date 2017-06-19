function plot_param_dis(obj,whichOne)
% plot the parameter disturbution

binNum = 100;
inputsNames = obj.inputsNames;
outputNames = obj.outputNames;

figure;

switch whichOne
    case 'inputs'
        inNum = length(inputsNames);
        for j=1:inNum
            temp = obj.get_param(whichOne{1,j});
            subplot(ceil(inNum/2),floor(inNum/2),j)
            hist(temp,binNum);    grid minor;
            title(['historiogram of ',whichOne{1,j}]);
        end
    case 'outputs'
        outNum = length(outputNames);
        for j=1:outNum
            temp = obj.get_param(whichOne{1,j});
            subplot(ceil(outNum/2),floor(outNum/2),j)
            hist(temp,binNum);    grid minor;
            title(['historiogram of ',whichOne{1,j}]);
        end
    otherwise
        temp = obj.get_param(whichOne);
        hist(temp,binNum);    grid minor;
        title(['historiogram of ',whichOne]);
end

end

