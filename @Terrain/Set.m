function [Te] = Set(Te,varargin)
% Sets desired object properties
% Use as Set(Te,'Type',1,'init_slope',0,_)

nParams = (nargin-1)/2;
if rem(nParams,1)~=0 || nargin<1
    error('Set failed: not enough inputs')
else
    for p = 1:nParams
        key = varargin{2*p-1};
        value = varargin{2*p};
        if ~isnumeric(value) && ~strcmp(key,'Type')
            error('Set failed: property value must be numeric');
        end

        switch key
            case 'Type'
                if isnumeric(value)
                    Te.Type = value;
                else
                    switch value
                        case {'inc','inclined','inclined plane'}
                            Te.Type = 0;
                        case {'sine','sinus','sinusoidal'}
                            Te.Type = 1;
                        case {'inf','infinite','inf par','infpar','infinite parabolla'}
                            Te.Type = 2;
                        case {'fin','finite','finite par','finite parabolla'}
                            Te.Type = 3;
                        otherwise
                            error('Set failed: Terrain type unknown');
                    end
                end
            case 'sinAmp'
                Te.sinAmp = value;
                Te.Type = 1;
            case 'sinFreq'
                Te.sinFreq = value;
                Te.Type = 1;
            case 'parK' % Parabolla constant parK/2*x^2
                Te.parK = value;
            case 'start_slope'
                Te.start_slope = value;
            case 'end_slope'
                Te.end_slope = value;
            case 'start_x'
                Te.start_x = value;
            case 'end_x'
                Te.end_x = value;
            case 'FloorStep'
                Te.FloorStep = value;
            case 'VertLines'
                Te.VertLines = value;
                Te.FloorVLx=zeros(1,NLines);
                Te.FloorVL=zeros(1,NLines);
            case 'FloorColor'
                Te.FloorColor = value;
            case 'LineWidth'
                Te.LineWidth = value;
        end
    end
    Te=SetEndConditions(Te);
end
