function [vec_shuffled] = shuffle_vec(vec)
% SHUFFLE_VEC takes a vector and shuffles its elements

vec_shuffled = vec(randperm(length(vec)));
end

