function [ NC ] = AddPulse( NC, varargin )
% Adds a pulse to the controller object
% Pulse parameters required:
% joint - Which joint to apply the pulse to
% amp - Pulse amplitude
% offset - Defines beginning of pulse as % of neuron period
% dur - Pulse duration as % of neuron period

% If the individual feedback is selected then the
% gains (k_u - up and k_d - down) can also be provided

nParams = (nargin-1)/2;
if rem(nParams,1)~=0 || nargin<1
    error('Set failed: not enough inputs')
else
    % Get ready to add a new pulse
    PulID = NC.nPulses + 1;
    
    for p = 1:nParams
        key = varargin{2*p-1};
        value = varargin{2*p};
        if ~isnumeric(value)
            error('Set failed: property value must be numeric');
        end
        
        switch key
            case 'joint'
                NC.OutM(value,PulID) = 1;    
            case 'amp'
                NC.Amp0(PulID) = value;
                NC.Amp(PulID) = value;
            case 'offset'
                NC.Offset(PulID) = value;
            case 'dur'
                NC.Duration(PulID) = value;
            case 'k_u'
                NC.kTorques_u(PulID) = value;
            case 'k_d'
                NC.kTorques_d(PulID) = value;

            otherwise
                error(['Set failed: ',key,' property not found']);
        end
    end
    
    % Update number of pulses
    NC.nPulses = PulID;
    NC.nEvents = 2 + 2*PulID;
    NC.Switch(PulID,1) = 0;

end

