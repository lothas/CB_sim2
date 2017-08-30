function fh = MoE_calc_fh(obj,inputs,targets,...
    expertsNN,g)
% this function run all of the data on the experts, then calculate the
% squar error. it also calculate the "fh" for the MoE training.

% TODO: find a way to make it work with cases 9 and 10 of our paper 
% % % % 2 parameters in the NN output.
num_of_train_samples = size(targets,2);

errMat = zeros(obj.expertCount,num_of_train_samples);
for j=1:obj.expertCount % run the data throught the experts to get initial clustering
    tempNet = expertsNN{1,j};
    outMat = tempNet(inputs);
    errMat(j,:) = outMat - targets;
end
seMat = errMat.^2; % squar error

% calc f_h from the equation in Jacobs1990 paper
% as seMat = yStar_yi in the paper.
fh = zeros(obj.expertCount,num_of_train_samples);
for k=1:num_of_train_samples
    fh(:,k) = g(:,k) .* exp(-0.5 .* seMat(:,k) );
    fh(:,k) = fh(:,k) ./ sum(fh(:,k),1);
end

% % % check if the sum of the gait for each sample is '1'
% if (sum(sum(fh,1),2) ~=  num_of_train_samples)
%    error('something is wrong with "f_h"!');
% end

end

