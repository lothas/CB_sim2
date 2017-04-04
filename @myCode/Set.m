function [obj] = Set( obj, varargin )
% Sets desired object properties

switch varargin{1}
    case {'NN','normalNN'}
        obj.NN.hiddenNeuronNum = varargin{2};
        obj.NN.net = feedforwardnet(varargin{2});
        obj.numOfIteretions = varargin{3};
        obj.NN.net.trainParam.epochs = varargin{3};
        
    case {'my_MoE','our_MoE'}
        if nargin == 6
            obj.numOfIteretions = varargin{2}; 
            obj.expertCount = varargin{3};
            obj.my_MoE_out.ExpertHidLayer = varargin{4};
            obj.my_MoE_out.maxEphocs = varargin{5};
            obj.my_MoE_out.perf_Stop_cond = 1e-7;
            obj.my_MoE_out.gradient_stop = 1e-5;
            
        else
           disp(['input variable order:',...
               'numOfIteretions,expertCount,ExpertHidLayer,GateHidLayer',...
               'maxEphocs,competetiveFlag']);
           error('not enougth inputs variable');
        end
        
    case {'paper_MoE','Jordan_papers'}
        obj.numOfIteretions = varargin{2};
        obj.expertCount = varargin{3};
        obj.paper_MoE_out.learningRate = varargin{4};
        obj.paper_MoE_out.decay = varargin{5};
    otherwise
        error('input not valid!');
end

% if we use MoE, create colors vewctor and legend string array:
if ~isempty(obj.expertCount) && (obj.expertCount > 1)
    % set the experts colors for the graph colors 
    switch obj.expertCount
        case {2,3} % in case of small number of expert, make colors clear:
            obj.colors = [1,0,0;0,1,0;0,0,1]; % RGB
        otherwise
            obj.colors = rand(obj.expertCount,3); % random colors
    end
    % defining the Legend with the experts names
    legendNames_temp = cell(1,obj.expertCount);
    for i=1:obj.expertCount
        legendNames_temp{1,i} = ['#',num2str(i),' expert'];
    end
    obj.legendNames = legendNames_temp;
end

end

