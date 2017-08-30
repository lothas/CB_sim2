function obj = MoE_train(obj)
% this function calls the different Mixute of Experts training functions

switch obj.MoE_method
    case 'collaboration'
        obj = MoE_train_collaboration(obj);
    case {'hardCompetetetive','softCompetetetive'}
        obj = MoE_train_Competetive(obj);
    otherwise
        error('invalid "MoE" Method');
        
end

end

