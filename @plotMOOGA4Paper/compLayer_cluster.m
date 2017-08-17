function [net,idx] = compLayer_cluster(X1,X2)
% clustering using competetive layer NN

% Inputs:
% *) 'X1' 'X2' - inputs to NN
% 
% Outputs:
% *) 'net' - NN object
% *) 'idx' - clustering indices
% 

inputs = ([X1,X2])';

net = competlayer(num_of_clusters);
net = train(net,inputs);
view(net)
outputs = net(inputs);
idx = vec2ind(outputs);

end

