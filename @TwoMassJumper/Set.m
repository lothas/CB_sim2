function [ TMJ ] = Set( TMJ, varargin )
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
            case 'm1' % upper mass
                TMJ.m1 = value;
            case 'm2' % lower mass
                TMJ.m2 = value;
            case 'k' % spring constant
                TMJ.k = value;
            case 'l0' % spring zero length
                TMJ.l0 = value;
            case 'g' % gravity
                TMJ.g = value;
            case 'damp' % Damping
                TMJ.damp = value;
            case 'Phase' % Phase (stance or flight)
                TMJ.Phase = value;
                
            % %%%%%% % Render parameters % %%%%%% %
            case 'm_size'
                TMJ.m_size = value;
            case 'm_color'
                TMJ.m_color = value;
            case 'LineWidth'
                TMJ.LineWidth = value;
            otherwise
                error(['Set failed: ',key,' property not found']);
        end
        
        TMJ = TMJ.SetND();
    end
end

