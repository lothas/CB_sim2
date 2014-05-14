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
    PulseParams = cell(6,1);
        
    for p = 1:nParams
        key = varargin{2*p-1};
        value = varargin{2*p};
        
        switch lower(key)
            case 'joint'
                PulseParams{1} = value;
            case 'amp'
                PulseParams{2} = value;
            case 'offset'
                PulseParams{3} = value;
            case 'dur'
                PulseParams{4} = value;
            case 'k_u'
                PulseParams{5} = value;
            case 'k_d'
                PulseParams{6} = value;

            otherwise
                error(['Set failed: ',key,' property not found']);
        end
    end
    
    % Check if all the required parameters where input
    if NC.FBType == 2
        Missing = find(cellfun('isempty',PulseParams)==1);
    else
        Missing = find(cellfun('isempty',PulseParams(1:4))==1);
    end
    
    if isempty(Missing)
        % Add the new pulse
        % Is it an externally triggered pulse?
        if ~isnumeric(PulseParams{3}) % should be 'ext'
            % Add external pulse
            NC.ExtPulses = [NC.ExtPulses, PulID];
            PulseParams{3} = 2; % will never be reached
        end
        
        if any(cellfun(@isnumeric,PulseParams([1,2,4:end]))==0)
            error('Set failed: property value must be numeric');
        end
        
        NC.OutM(PulseParams{1},PulID) = 1;
        NC.Amp0(PulID) = PulseParams{2};
        NC.Amp(PulID) = PulseParams{2};
        NC.Offset(PulID) = PulseParams{3};
        NC.Duration(PulID) = PulseParams{4};
        if NC.FBType == 2
            NC.kTorques_u(PulID) = PulseParams{5};
            NC.kTorques_d(PulID) = PulseParams{6};
        end

        % Update number of pulses
        NC.nPulses = PulID;
        NC.nEvents = 2 + 2*PulID;
        NC.Switch(PulID,1) = 0;
    else
        Keys = {'Joint', 'Amp', 'Offset', 'Dur', 'k_u', 'k_d'};
        error(['No value provided for ',Keys{Missing},' parameter(s)']);
    end
end

