classdef Genome
    %GENOME Holds the keys to encoding/decoding a genome
    %   An instance of a Genome object uses a set of keys
    %   to decode a string of parameters and set the
    %   simulation properties (Model, Controller, etc)
    
    properties
        Keys
        KeyLength
        KeyExtra
        Length
    end
    
    methods
        function Ge = Genome(varargin)
            switch nargin
                case 1
                    Ge = Ge.SetKeys(varargin{1});
            end
            
            % KeyLengths can be set here or from the outside
            % They are used to define the number of genes for a
            % specific key
            Ge.KeyLength.pulse = 3;
            
            % KeyExtras can be set here or from the outside
            % They are used to provide extra information about
            % how to decode a key, e.g. where to put every
            % initial condition given
            Ge.KeyExtra.IC = [-1,1,2,3,4];
            % Here ^ 4 values should be provided for an IC vector
            % with 5 coordinates (the first value is substractred
            % from the first IC)
        end
        
        function Ge = SetKeys(Ge, Keys)
            Ge.Keys = Keys;
            Ge.Length = 0;
            for k = 1:size(Keys, 2)
                if isfield(Ge.KeyLength,Keys{1,k})
                    Ge.Length = Ge.Length + ...
                        Ge.KeyLength.(Keys{1,k})*Keys{2,k};
                else
                    Ge.Length = Ge.Length + Keys{2,k};
                end
            end
        end
    end
    
end

