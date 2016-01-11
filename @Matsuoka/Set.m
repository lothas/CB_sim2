function [ MO ] = Set( MO, varargin )
% Sets desired object properties
% Use as Set(CB,'m',3,'mh',10,_)

nParams = (nargin-1)/2;
if rem(nParams,1)~=0 || nargin<1
    error('Set failed: not enough inputs')
else
    for p = 1:nParams
        key = varargin{2*p-1};
        value = varargin{2*p};
        if ~isnumeric(value)
            error('Set failed: property value must be numeric');
        end
        
        switch lower(key)
            case {'npulses','nneurons','n_pulses','n_neurons'}
                MO.nPulses = value; % Overall number of E/F neuron pairs
                
            % neuron parameters
            case {'tau', 'tau_u'}
                MO.tau = value;
            case {'tav', 'tau_v'}
                MO.tav = value;
            case 'beta'
                MO.beta = value;
            case 'win' % Intra-neuron pair inhibition
                N = length(value);
                if N == 1
                    % One inhibition to rule them all
                    D = mat2cell(repmat([0, value; value, 0], ...
                                 1, MO.nPulses), 2, 2*ones(1,MO.nPulses)); 
                else
                    D = {};
                    if N == MO.nPulses
                        % One inhibition per neuron pair
                        for i = 1:length(value)
                            D = [D, [0, value(i); value(i), 0]];  %#ok<AGROW>
                        end
                    else
                        % One inhibition per neuron
                        for i = 1:length(value)
                            D = [D, [0, value(2*i-1); value(2*i), 0]];  %#ok<AGROW>
                        end
                    end
                end
                MO.win = blkdiag(D{:});
            case 'wex' % Extra-neuron inhibition (E to E and F to F)
                MO.wex = zeros(length(value));
                v = 1;
                for n = 1:length(value)
                    ids = 2-mod(n,2):2:2*MO.nPulses;
                    ids(ids == n) = [];
                    MO.wex(n,ids) = value(v:v+length(ids)-1);
                    v = v+length(ids);
                end
                
            % Controller Output
            case {'amp0', 'amp'} % Base neuron amplitude multiplier
                N = length(value);
                if N == MO.nPulses
                    MO.Amp0 = reshape([value; value], 1, []);
                else
                    MO.Amp0 = value;
                end
                MO.Amp = MO.Amp0;
                
            % Feedback
            case {'fbtype', 'feedback', 'fb'}
                MO.FBType = value;
                % Not yet implemented
                
            % Gains
            case {'ks_tau', 'speed_tau', 'tau_speed_gain'}
                MO.ks_tau = value;
            case {'ks_out', 'speed_out', 'torque_speed_gain'}
                MO.ks_out = value;
            
            otherwise
                error(['Set failed: ',key,' property not found']);
        end
    end
end

