function [outputs] = apply_net(obj,inputs,method)
% apply NN/MoE to get the outputs

switch method
    case 'NN'
        % train NN:
        net = obj.NN.net;
        outputs = net(inputs);
    case 'MoE colaboration'
        [outputs,~,~,~] =...
            obj.MoE.MoE_testNet(inputs,...
            obj.MoE.expertsNN,obj.MoE.gateNet);
    case 'MoE hard'
        [outputs,~,~,~] =...
            obj.MoE.MoE_testNet(inputs,...
            obj.MoE.expertsNN,obj.MoE.gateNet);
    case 'MoE soft'
        [outputs,~,~,~] =...
            obj.MoE.MoE_testNet(inputs,...
            obj.MoE.expertsNN,obj.MoE.gateNet);
    otherwise
        error('invalid method...');
end

end

