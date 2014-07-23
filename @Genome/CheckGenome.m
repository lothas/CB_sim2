function [ Result, Seq ] = CheckGenome( Ge, Seq )
% Checks a gene sequence to verify that it satisfies
% all the min/max conditions
    Result = {1,''};
    
    % Check that it's the correct length
    if length(Seq)~=Ge.Length
        Result = {0;
            ['Decode failed. Genes provided: ',num2str(length(Seq)),...
            '. Genes required: ',num2str(Ge.Length)]};
        return;
    end
    
    % Check that it's between the max and min range
    if any(Seq<Ge.Range(1,:) | Seq>Ge.Range(2,:))
        Result = {0;
        	'Sequence outside allowed genome range'};
        return;
    end

    SeqPos = 1;
    for k = 1:size(Ge.Keys,2)
        % Go over each key
        switch Ge.Keys{1,k}
            case 'Pulses'
                NPulses = Ge.Keys{2,k}(1);
                SeqPos0 = SeqPos;
                x = 0:0.001:1;
                Torque = zeros(1,length(x));
                for p = 1:NPulses
                    % Check that the pulse ends before the phase resets
                    End = Seq(SeqPos0+1)+Seq(SeqPos0+2);
                    if End>1
                        % Shorten both the offset and the duration
                        % so that End will be slightly smaller than 1
                        Seq(SeqPos0+1) = Seq(SeqPos0+1)/(End*1.001);
                        Seq(SeqPos0+2) = Seq(SeqPos0+2)/(End*1.001);
                        Result{2} = 'Pulse shortened';
                    end
                    
                    Torque = Torque + ...
                        Seq(SeqPos0)*(x>Seq(SeqPos0+1) & x<End);
                    SeqPos0 = SeqPos0 + Ge.KeyLength.Pulses;
                end
                
                % Check that the compound torque signal doesn't exeed
                % the min/max allowed
                factor = [];
                if any(Torque<Ge.Range(1,SeqPos)) % Min
                    factor = min(Torque)/Ge.Range(1,SeqPos);
                end
                if any(Torque>Ge.Range(2,SeqPos)) % Max
                    factor = max(Torque)/Ge.Range(2,SeqPos);
                end
                if ~isempty(factor)
                    % Fix the genome
                    SeqPos0 = SeqPos;
                    for p = 1:NPulses
                        Seq(SeqPos0) = Seq(SeqPos0)/factor;
                        SeqPos0 = SeqPos0 + Ge.KeyLength.Pulses;
                    end
                    Result{2} = ['Compound torque reduced by factor of ',num2str(factor)];
                end
            case 'ExtPulses'
                NPulses = Ge.Keys{2,k}(1);
                SeqPos0 = SeqPos;
                
                Total = 0;
                for p = 1:NPulses
                    Total = Total + Seq(SeqPos0);
                    SeqPos0 = SeqPos0 + Ge.KeyLength.ExtPulses;
                end
                
                factor = [];
                if Total<Ge.Range(1,SeqPos) % Min
                    factor = Total/Ge.Range(1,SeqPos);
                end
                if Total>Ge.Range(2,SeqPos) % Max
                    factor = Total/Ge.Range(2,SeqPos);
                end
                if ~isempty(factor)
                    % Fix the genome
                    SeqPos0 = SeqPos;
                    for p = 1:NPulses
                        Seq(SeqPos0) = Seq(SeqPos0)/factor;
                        SeqPos0 = SeqPos0 + Ge.KeyLength.ExtPulses;
                    end
                    Result{2} = ['Compound event triggered torque reduced',...
                        'by factor of ',num2str(factor)];
                end
        end
        
        % Move the sequence reading position
        SeqPos = Ge.AdvSeq(SeqPos,k);
    end
end

