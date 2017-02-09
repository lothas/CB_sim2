function [obj] = Set( obj, varargin )
% Sets desired object properties

switch varargin{1}
    case {'NN','normalNN'}
        obj.NN.hiddenNeuronNum = varargin{2};
        obj.NN.net = feedforwardnet(varargin{2});
    case {'my_MoE','our_MoE'}
        %
    case {'paper_MoE','Jordan_papers'}
        %
    otherwise
        error('input not valid!');
end

end

