function [samplNew] = matsuoka_arrange_uniq(obj,sampl)
% this function takes a Matsuoka CPG and arranges all of the parameters
% in a unique way based on the largest weight for every matsuoka neuron

% outputs:
% 1) 'samplNewUniq' - the inputs in a uniq way
% 2) 'TargetNewUniq' - the targets

numOfSAmple = size(sampl,2);
samplNew = zeros(size(sampl));

orderOfSampl = obj.inputsNames;
location_tau = strcmp('tau',orderOfSampl);
location_b = strcmp('b',orderOfSampl);
location_c1 = strcmp('c_1',orderOfSampl);
location_c2 = strcmp('c_2',orderOfSampl);
location_c3 = strcmp('c_3',orderOfSampl);
location_c4 = strcmp('c_4',orderOfSampl);
location_w12 = strcmp('w_{12}',orderOfSampl);
location_w13 = strcmp('w_{13}',orderOfSampl);
location_w14 = strcmp('w_{14}',orderOfSampl);
location_w21 = strcmp('w_{21}',orderOfSampl);
location_w23 = strcmp('w_{23}',orderOfSampl);
location_w24 = strcmp('w_{24}',orderOfSampl);
location_w31 = strcmp('w_{31}',orderOfSampl);
location_w32 = strcmp('w_{32}',orderOfSampl);
location_w34 = strcmp('w_{34}',orderOfSampl);
location_w41 = strcmp('w_{41}',orderOfSampl);
location_w42 = strcmp('w_{42}',orderOfSampl);
location_w43 = strcmp('w_{43}',orderOfSampl);

tau = sampl(location_tau,:);
b = sampl(location_b,:);
c_1 = sampl(location_c1,:);
c_2 = sampl(location_c2,:);
c_3 = sampl(location_c3,:);
c_4 = sampl(location_c4,:);
w_12 = sampl(location_w12,:);
w_13 = sampl(location_w13,:);
w_14 = sampl(location_w14,:);
w_21 = sampl(location_w21,:);
w_23 = sampl(location_w23,:);
w_24 = sampl(location_w24,:);
w_31 = sampl(location_w31,:);
w_32 = sampl(location_w32,:);
w_34 = sampl(location_w34,:);
w_41 = sampl(location_w41,:);
w_42 = sampl(location_w42,:);
w_43 = sampl(location_w43,:);

if (obj.sizeOfCPG == 2) && (isempty(w_12) || isempty(w_21))
    if obj.disp_information
        disp('not taining on all of the CPG weights, continuing without rearenging!');
    end
    samplNew = sampl;
    return
end

W_cond = (isempty(w_12) || isempty(w_13) || isempty(w_14) ||...
    isempty(w_21) || isempty(w_23) || isempty(w_24) ||...
    isempty(w_31) || isempty(w_32) || isempty(w_34) ||...
    isempty(w_41) || isempty(w_42) || isempty(w_43));
if (obj.sizeOfCPG == 4) && W_cond
    if obj.disp_information
        disp('not taining on all of the CPG weights, continuing without rearenging!');
    end
    samplNew = sampl;
    return
end

if ~isempty(tau)
    samplNew(location_tau,:) = tau;
end

if ~isempty(b)
    samplNew(location_b,:) = b;
end


switch obj.sizeOfCPG
    case 2 % if we have CPG with 2 neurons
        newSampl_Weights = zeros(2,numOfSAmple);
        newSampl_c = zeros(2,numOfSAmple);
        for i=1:numOfSAmple
            mat = [0        ,w_12(1,i);
                   w_21(1,i),0        ];

            biggestWeights = max(mat,[],1);
            [~,sortedWeightsIND] = sort(biggestWeights,2);

            mat = mat(sortedWeightsIND',sortedWeightsIND');

            newSampl_Weights(:,i) = [mat(1,2);mat(2,1)];
            if ~isempty(c_1) || ~isempty(c_2)                  
                C = [c_1(1,i);c_2(1,i)];
                newSampl_c(1,i) = C(sortedWeightsIND',1);
            end
        end
        locationOfWeights = location_w12 | location_w21;
        samplNew(locationOfWeights,:) = newSampl_Weights;
        if ~isempty(c_1) || ~isempty(c_2)
            locationOfCs= location_c1 | location_c2;
            samplNew(locationOfCs,:) = newSampl_c;
        end
        
    case 4 % if we have CPG with 4 neurons
        newSampl_Weights = zeros(12,numOfSAmple);
        newSampl_c = zeros(4,numOfSAmple);
        for i=1:numOfSAmple
            mat = [0        ,w_12(1,i),w_13(1,i),w_14(1,i);
                   w_21(1,i),0        ,w_23(1,i),w_24(1,i);
                   w_31(1,i),w_32(1,i),0        ,w_34(1,i);
                   w_41(1,i),w_42(1,i),w_43(1,i),0       ];

            biggestWeights = max(mat,[],1);
            [~,sortedWeightsIND] = sort(biggestWeights,2);

            mat = mat(sortedWeightsIND',sortedWeightsIND');

            newSampl_Weights(:,i) = [mat(1,2);mat(1,3);mat(1,4);...
                                mat(2,1);mat(2,3);mat(2,4);...
                                mat(3,1);mat(3,2);mat(3,4);...
                                mat(4,1);mat(4,2);mat(4,3)];
            if ~isempty(c_1) || ~isempty(c_2) || ~isempty(c_3) || ~isempty(c_4)                  
                C = [c_1(1,i);c_2(1,i);c_3(1,i);c_4(1,i)];
                newSampl_c(1,i) = C(sortedWeightsIND',1);
            end
        end
        locationOfWeights = location_w12 | location_w13 | location_w14 |...
            location_w21 | location_w23 | location_w24 |...
            location_w31 | location_w32 | location_w34 |...
            location_w41 | location_w42 | location_w43;
        samplNew(locationOfWeights,:) = newSampl_Weights;
        if ~isempty(c_1) || ~isempty(c_2) || ~isempty(c_3) || ~isempty(c_4)
            locationOfCs= location_c1 | location_c2 | location_c3 | location_c4;
            samplNew(locationOfCs,:) = newSampl_c;
        end
        
    otherwise
        error('size of CPG is unknown!');
end

end

