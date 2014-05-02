function [ NC ] = Set( NC, varargin )
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
        
        switch key
            % LIF parameters
            case 'P_reset'
                NC.P_reset = value;
            case 'P_th'
                NC.P_th = value;
            case 'P_0' % Initial phase for oscillator
                NC.P_0 = value;
            case 'P_LegE'
                NC.P_LegE = value; % timed leg extension
            case 'omega0' % Oscillator's frequency
                NC.omega0 = value;
                NC.omega = value;
                
            % Controller Output
            case 'nPulses'
                NC.nPulses = value;  % Overall number of pulses
            case 'Amp0'
                NC.Amp0 = value;     % Base pulse amplitude in N
                NC.Amp = value;      % Pulse amplitude in N
            case 'Offset'
                NC.Offset = value;   % Defines beginning of pulse as % of neuron period
            case 'Duration'
                NC.Duration = value; % Pulse duration as % of neuron period
            case 'pSoff'
                NC.pSoff = value;    % Phase at which to turn off external inputs

            % Adaptation
            case 'FBType'
                NC.FBType = value;
                % 0 - no feedback, 1 - single gain for omega and each joint
                % 2 - individual gains for each pulse

            % Gains
            case 'kOmega_u'
                NC.kOmega_u = value;
            case 'kOmega_d'
                NC.kOmega_d = value;
            case 'kTorques_u'
                NC.kTorques_u = value;
            case 'kTorques_d'
                NC.kTorques_d = value;
            
            otherwise
                error(['Set failed: ',key,' property not found']);
        end
    end
end

