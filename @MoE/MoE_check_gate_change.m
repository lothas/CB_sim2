function [max_err_g,max_err_g_ind,err_g] = MoE_check_gate_change(obj,g_old,g_new)
% this function is called every training iteration and is checking, for
% every point, the change in the gate output.

% OUTPUTS:
%	*) max_err_g - the maximum gate change
%	*) max_err_g_ind - the index to the sample with the maximum gate change
%	*) err_g - a row vector contain the change for all of the samples

% calculate change by mse error:
err_g = zeros(1,size(g_old,2));
for j=1:size(g_old,2)
    err_g(1,j) = immse(g_old(:,j),g_new(:,j));
end
%    check for which point we get the maximum change.
[max_err_g,max_err_g_ind] = max(err_g,[],1);

% % option 2: calculate change by calculating the norm of the difference
% % between old to new.
% delta_g = g_old - g_new;
% err_g2 = zeros(1,size(delta_g,2));
% for j=1:size(delta_g,2)
%     err_g2(1,j) = norm(delta_g(:,j));
% end

end

