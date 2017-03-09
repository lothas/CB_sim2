function [obj] = Set( obj, varargin )
% Sets desired object properties

switch varargin{1}
    case {'NN','normalNN'}
        obj.NN.hiddenNeuronNum = varargin{2};
        obj.NN.net = feedforwardnet(varargin{2});
        obj.numOfIteretions = varargin{3};
        obj.NN.net.trainParam.epochs = varargin{3};
        
    case {'my_MoE','our_MoE'}
        if nargin == 8
            obj.numOfIteretions = varargin{2};
            obj.expertCount = varargin{3};
            obj.my_MoE_out.ExpertHidLayer = varargin{4};
            obj.my_MoE_out.GateHidLayer = varargin{5};
            obj.my_MoE_out.maxEphocs = varargin{6};
            obj.my_MoE_out.competetiveFlag = varargin{7};
            % set the experts colors for the graph colors 
            switch varargin{3} %obj.expertCount
                case {2,3} % in case of small number of expert, make colors clear:
                    obj.colors = [1,0,0;0,1,0;0,0,1]; % RGB
                otherwise
                    obj.colors = rand(varargin{3},3); % random colors
            end
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
        % set the experts colors for the graph colors 
        switch varargin{3} %obj.expertCount
            case {2,3} % in case of small number of expert, make colors clear:
                obj.colors = [1,0,0;0,1,0;0,0,1];
            otherwise
                obj.colors = rand(varargin{3},3);
        end
        
    otherwise
        error('input not valid!');
end

end

