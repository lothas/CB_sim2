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
                MO.tau0 = value;
            case {'tav', 'tau_v'}
                MO.tav = value;
                MO.tav0 = value;
            case {'tau_ratio','\tau_ratio'}
                MO.tau_ratio = value;
            case {'tau_r','\tau_r'}
                MO.tau = value;
                MO.tau0 = value;
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
            case '2neuron_symm_weights' % Neuron connection weights, symmetric 2neuron CPG
                MO.win = 0;
                MO.wex = zeros(2*MO.nPulses);
                % MO.wex = [0,value;value,0];
                % only the hip neurons are connects to each other
                MO.wex = [0       ,0       ,0       ,0       ;
                          0       ,0       ,0       ,0       ;
                          0       ,0       ,0       ,value   ;
                          0       ,0       ,value   ,0       ];
            case '2neuron_general_weights' % Neuron connection weights, general 2neuron CPG
                MO.win = 0;
                MO.wex = zeros(2*MO.nPulses);
                % MO.wex = [0,value(1);value(2),0];
                MO.wex = [0       ,0       ,0       ,0       ;
                          0       ,0       ,0       ,0       ;
                          0       ,0       ,0       ,value(1);
                          0       ,0       ,value(2),0       ];
            case '4neuron_symm_weights' % Neuron connection weights, symmetric 4neuron CPG
                MO.win = 0;
                MO.wex = zeros(2*MO.nPulses);
                MO.wex = [0       ,value(1),value(2),value(3);
                          value(1),0       ,value(4),value(5);
                          value(2),value(4),0       ,value(6);
                          value(3),value(5),value(6),0       ];
            case '4neuron_taga_like'
                MO.win = 0;
                MO.wex = zeros(2*MO.nPulses);
                MO.wex = [0       ,value(1),0       ,0;
                          value(1),0       ,value(2),value(3);
                          0       ,0       ,0       ,value(4);
                          0       ,0       ,value(4),0       ];
                      
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
                
            case 'amp_2n_same_inputs' % same tonic inputs to both MN
                MO.Amp0 = [0;0;value(1);value(1)];
                MO.Amp = MO.Amp0;
            case 'amp_2n_dif_inputs' % different tonic inputs to both MN
                MO.Amp0 = [0;0;value(1);value(2)];
                MO.Amp = MO.Amp0;
                
            case 'amp_4n_symm'
                MO.Amp0 = [value;value;value;value];
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
            case {'ks_c_2n_symm'} % when we have only the hip joint
                MO.ks_out = [0;0;value;value];
            case {'ks_c_2n_general'} % when we have only the hip joint
                MO.ks_out = [0;0;value(1);value(2)];
            case {'ks_c_4n_symm'}
                MO.ks_out = [value;value;value;value];
            
            otherwise
                error(['Set failed: ',key,' property not found']);
        end
    end
end

