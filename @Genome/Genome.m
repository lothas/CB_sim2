classdef Genome
    %GENOME Holds the keys to encoding/decoding a genome
    %   An instance of a Genome object uses a set of keys
    %   to decode a string of parameters and set the
    %   simulation properties (Model, Controller, etc)
    
    properties
        % Sequence decoding parameters
        Keys
        KeyLength
        Segments
        KeyExtra
        Length
        
        % Sequence range
        Range % Min/Max range for each gene
        
        % Crossover parameters
        COType = '2point';
        
        % Mutation parameters
        MutProb = 0.5;      % Probability that a single gene will mutate
        MutDelta = 0.1;     % Max strength of mutation as percentage of range
        MutType = 'uni';    % Type of mutation: uniform or normal
    end
    
    methods
        function Ge = Genome(varargin)            
            % KeyLengths can be set here or from the outside
            % They are used to define the number of genes for a
            % specific key
            Ge.KeyLength.Pulses = 3;
            Ge.KeyLength.ExtPulses = 2;
            
            % KeyExtras can be set here or from the outside
            % They are used to provide extra information about
            % how to decode a key, e.g. where to put every
            % initial condition given
            Ge.KeyExtra.IC = [1,-1,2,3,4];
            % Here ^ 4 values should be provided for an IC vector
            % with 5 coordinates (the first value is substractred
            % from the first IC)
            
            switch nargin
                case 1
                    Ge = Ge.SetKeys(varargin{1});
                case 2
                    Ge = Ge.SetKeys(varargin{1});
                    Ge = Ge.SetRange(varargin{2});
                case 3
                    Ge.KeyLength = varargin{2};
                    Ge = Ge.SetKeys(varargin{1});
                    Ge = Ge.SetRange(varargin{3});
            end
        end
        
        function Ge = SetKeys(Ge, Keys)
            Ge.Keys = Keys;
            Ge.Length = 0;
            Ge = Ge.SetSegments();
            Ge.Length = sum(Ge.Segments);
        end
        
        function Ge = SetRange(Ge, Range)
            Ge.Range = zeros(2,Ge.Length);
            SeqPos = 1; % Position along the genome sequence
            for r = 1:size(Range, 2)
                thisRange = cell2mat(Range(:,r));
                if isfield(Ge.KeyLength,Ge.Keys{1,r})
                    k = Ge.Keys{2,r}(1);
                    thisRange = repmat(thisRange,1,k);
                end
                tRL = size(thisRange,2);
                Ge.Range(:,SeqPos:SeqPos+tRL-1) = thisRange;
                SeqPos = SeqPos+tRL;
            end
        end
        
        function Ge = SetSegments(Ge)
            NKeys = size(Ge.Keys, 2);
            Ge.Segments = zeros(1,NKeys);
            for k = 1:NKeys
                if isfield(Ge.KeyLength,Ge.Keys{1,k})
                    Ge.Segments(k) = ...
                        Ge.KeyLength.(Ge.Keys{1,k})*Ge.Keys{2,k}(1);
                else
                    Ge.Segments(k) = Ge.Keys{2,k}(1);
                end
            end
        end
        
        function SeqPos = AdvSeq(Ge,SeqPos,k)
            SeqPos = SeqPos + Ge.Segments(k);
        end
        
        function Seq = SwitchDir(Ge,Seq)
            SeqPos = 1;
            for k = 1:size(Ge.Keys,2)
                % Go over each key
                switch Ge.Keys{1,k}
                    case {'Pulses','ExtPulses'}
                        for g = 1:Ge.Keys{2,k}(1)
                            Gene = SeqPos + ...
                                (g-1)*Ge.KeyLength.(Ge.Keys{1,k});
                            Seq(Gene) = -Seq(Gene);
                        end
                    case {'AngVelImp'}
                        Seq(SeqPos:SeqPos+Ge.Segments(k)-1) = ...
                            -Seq(SeqPos:SeqPos+Ge.Segments(k)-1);
                    case {'ks_out','amp'}
                        temp = Seq(SeqPos:SeqPos+Ge.Segments(k)-1);
                        for g = 1:Ge.Segments(k)/2
                            % Switch flexor and extensor signal amplitude
                            Gene = SeqPos + (g-1)*2;
                            Seq(Gene) = Seq(Gene+1);
                            Seq(Gene+1) = temp(1+(g-1)*2);
                        end 
                end
                SeqPos = SeqPos + Ge.Segments(k);
            end
        end
        
        function [seqStr] = seq2str(Ge, seq, precision) %#ok<INUSL>
            if nargin<3
                precision = 6;
            end
            cut = floor(50/precision);
            seqStr = ['seq = [',num2str(seq(1),precision)];
            for g = 2:length(seq)
                if mod(g,cut) == 0
                    seqStr = [seqStr,',...',10,num2str(seq(g),precision)]; %#ok<AGROW>
                else
                    seqStr = [seqStr,', ',num2str(seq(g),precision)]; %#ok<AGROW>
                end
            end
            seqStr = [seqStr,'];'];
        end
        
        function [genes] = GetGenes(Ge, seq, sel_keys)
            genes = seq;
            ng = 1;
            for i = 1:numel(sel_keys)
                key_id = find(strcmp(sel_keys{i},Ge.Keys(1,:)));
                if isempty(key_id)
                    error(['Key "',sel_keys{i},'" not found!']);
                end
                start_id = sum(Ge.Segments(1:key_id-1))+1;
                genes(ng:ng+Ge.Segments(key_id)-1) = ...
                    seq(start_id:start_id+Ge.Segments(key_id)-1);
                ng = ng+Ge.Segments(key_id);
            end
            genes = genes(1:ng-1);
        end
    end
    
end

