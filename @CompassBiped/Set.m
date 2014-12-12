function [ CB ] = Set( CB, varargin )
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
            case 'm' % leg mass
                k = CB.m_radius/CB.m;
                CB.m = value;
                CB.m_radius=k*CB.m;
            case 'mh' % hip mass
                k = CB.mh_radius/CB.mh;
                CB.mh = value;
                CB.mh_radius=k*CB.mh;
            case 'L' % leg length
                CB.L = value;
            case 'a' % leg center of mass
                CB.a = value;
            case 'I' % leg moment of inertia
                CB.I = value;
            case 'Clearance'
                CB.Clearance = value;
            case 'A2T'
                CB.A2T = value;
            case 'A2H'
                CB.A2H = value;
            case 'g' % gravity
                CB.g = value;
            case 'damp' % Damping (for both joints)
                CB.dampH = value;
                CB.dampA = value;
            case 'dampH' % Damping (for hip joint)
                CB.dampH = value;
            case 'dampA' % Damping (for ankle joint)
                CB.dampA = value;
            case 'xS' % Stance leg x position
                CB.xS = value;
            case 'yS' % Stance leg y position
                CB.yS = value;
            case 'Support' % Support leg (Left or Right)
                CB.Support = value;
                
            % %%%%%% % Render parameters % %%%%%% %
            case 'm_radius'
                CB.m_radius = value*CB.m;
            case 'm_color'
                CB.m_color = value;
            case 'mh_radius'
                CB.mh_radius = value*CB.mh;
            case 'mh_color'
                CB.mh_color = value;
            case 'leg_width'
                CB.leg_width = value;
            case 'leg_color'
                CB.leg_color = value;
            case 'CircRes'
                CB.CircRes = value;
            case 'LinkRes'
                CB.LinkRes = value;
            case 'LineWidth'
                CB.LineWidth = value;
            otherwise
                error(['Set failed: ',key,' property not found']);
        end
        
        CB = CB.SetND();
    end
end

