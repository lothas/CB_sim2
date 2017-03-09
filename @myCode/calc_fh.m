function [fh] = calc_fh(obj,expertsNN,g)
% this function run all of the data on the experts, then calculate the
% squar error. it also calculate the "fh" for the MoE training.


num_of_train_samples = size(obj.sampl_train,2);

errMat = zeros(obj.expertCount,num_of_train_samples);
for j=1:obj.expertCount % run the data throught the experts to get initial clustering
    tempNet = expertsNN{1,j};
    outMat = tempNet(obj.sampl_train);
    errMat(j,:) = outMat - obj.targ_train;
end
seMat = errMat.^2; % squar error

% calc f_h from the equation in Jacobs1990 paper
% as seMat = yStar_yi in the paper.
fh = zeros(obj.expertCount,num_of_train_samples);
for k=1:num_of_train_samples
    fh(:,k) = g(:,k) .* exp(-0.5 .* seMat(:,k) );
    fh(:,k) = fh(:,k) ./ sum(fh(:,k),1);
end

end

