function prepareData_NN(obj)
% prepare Matrices to NN training
% each ROW = different param.
% each Col = different sample.

for k=1:2
    if k == 1 %define inputs
        tempNames = obj.inputsNames;
    else %define outputs
        tempNames = obj.outputNames;
    end
    
    temp = zeros(length(tempNames),length(obj.ids));
    for i=1:length(tempNames)
        temp(i,:) = obj.get_param(tempNames{1,i});
    end
    
    if k == 1 %define inputs
        sampl = temp;
        clear temp
    else %define outputs
        targ = temp;
    end
end

obj.sampl = sampl;
obj.targ = targ;
end

