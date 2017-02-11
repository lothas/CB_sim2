function [obj] = Set( obj, varargin )
% Sets desired object properties

switch varargin{1}
    case {'NN','normalNN'}
        obj.NN.hiddenNeuronNum = varargin{2};
        obj.NN.net = feedforwardnet(varargin{2});
        
    case {'my_MoE','our_MoE'}
        if nargin == 8
            obj.numOfIteretions = varargin{2};
            obj.expertCount = varargin{3};
            obj.my_MoE_out.ExpertHidLayer = varargin{4};
            obj.my_MoE_out.GateHidLayer = varargin{5};
            obj.my_MoE_out.maxEphocs = varargin{6};
            obj.my_MoE_out.competetiveFlag = varargin{7};
        else
           disp(['input variable order:',...
               'numOfIteretions,expertCount,ExpertHidLayer,GateHidLayer',...
               'maxEphocs,competetiveFlag']);
           error('not enougth inputs variable');
        end
        
    case {'paper_MoE','Jordan_papers'}
        %
    otherwise
        error('input not valid!');
end

end

