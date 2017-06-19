function [normalize_data,normParams] = normData(obj,data)

    % normalize the data by ( x_norm = (x-x_mean)/stdev(x) )  

    normParams = zeros(size(data, 1), 2);
    normalize_data = zeros(size(data));
    for i = 1:size(data, 1)
        feat = data(i, :);
        normParams(i, :) = [mean(feat), std(feat)];
        normalize_data(i, :) = (feat - normParams(i, 1))/normParams(i, 2);
    end
end


