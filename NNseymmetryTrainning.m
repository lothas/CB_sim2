

%% cheking symmetry with uniqeness

    
newSampl = zeros(size(samplNewSymmetric,1),length(TargetNewSymmetric));

for i=1:length(TargetNewSymmetric)

    newSampl(1:2,i) = samplNewSymmetric(1:2,i);

    mat = [0             ,samplNewSymmetric(7,i)  ,samplNewSymmetric(8,i) ,samplNewSymmetric(9,i);
        samplNewSymmetric(10,i),0                 ,samplNewSymmetric(11,i),samplNewSymmetric(12,i);
        samplNewSymmetric(13,i),samplNewSymmetric(14,i) ,0                ,samplNewSymmetric(15,i);
        samplNewSymmetric(16,i),samplNewSymmetric(17,i) ,samplNewSymmetric(18,i),0               ];
    
    biggestWeights = max(mat,[],1);
    [~,sortedWeightsIND] = sort(biggestWeights,2);
    newSampl(3:6,i) = samplNewSymmetric(sortedWeightsIND+2,i);
    
    mat = mat(sortedWeightsIND',sortedWeightsIND');

    newSampl(7:18,i) = [mat(1,2);mat(1,3);mat(1,4);...
                        mat(2,1);mat(2,3);mat(2,4);...
                        mat(3,1);mat(3,2);mat(3,4);...
                        mat(4,1);mat(4,2);mat(4,3)];


end
% 
firstPerm = newSampl(1,1:dataPointsNum);
for i=1:24
    check = sum(newSampl(1,(dataPointsNum*i+1):dataPointsNum*(i+1)) == firstPerm);
    if check ~= dataPointsNum
        disp(['check fail: permutation ' , num2str(i), 'failed']);
    else
        disp(['check ',num2str(i), ' passed'])
    end
end