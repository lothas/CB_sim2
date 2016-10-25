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
            case {'tau_ratio','\tau_ratio'}
                MO.tau_ratio = value;
            case {'tau_r','\tau_r'}
                MO.tau = value;
                MO.tav = MO.tau_ratio*value;
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
                MO.wex = zeros(2*MO.nPulses);
                v = 1;
                for n = 1:2*MO.nPulses
                    ids = 1:2*MO.nPulses;
                    i0 = floor((n-1)/2)*2+1;
                    ids(i0:i0+1) = [];
                    MO.wex(n,ids) = value(v:v+length(ids)-1);
                    v = v+length(ids);
                end
            case 'weights' % Neuron connection weights, general
                MO.win = 0;
                MO.wex = zeros(2*MO.nPulses);
                v = 1;
                for i = 1:2*MO.nPulses
                    % Genes affect the coupling weight from neuron j to
                    % neuron i, for all j~=i
                    ids = 1:2*MO.nPulses;
                    ids(i) = [];
                    
                    MO.wex(i,ids) = value(v:v+length(ids)-1);
                    v = v+length(ids);
                end
                
            % Controller Output
            case {'amp0', 'amp', 'c_i'} % Base neuron amplitude multiplier
                N = length(value);
                if N == MO.nPulses
                    MO.Amp0 = reshape([value; value], [], 1);
%                     MO.Amp0 = reshape([value; value], 1, []);
                else
                    MO.Amp0 = value';
%                     MO.Amp0 = value;
                end
                MO.Amp = MO.Amp0;
                
            % Feedback
            case {'fbtype', 'feedback', 'fb'}
                MO.FBType = value;
                % Not yet implemented
                
            % Gains
            case {'ks_tau', 'speed_tau', 'tau_speed_gain', 'ks_\tau'}
                MO.ks_tau = value;
            case {'ks_out', 'speed_out', 'torque_speed_gain', 'ks_c'}
                MO.ks_out = value';
            
            otherwise
                error(['Set failed: ',key,' property not found']);
        end
    end
end

