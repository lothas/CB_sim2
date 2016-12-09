function dist = GetDistance( GA, Gen, type )
%GETDISTANCE Summary of this function goes here
%   Detailed explanation goes here
    if nargin<3
        type = 'seq';
    end
    
%     % Get all sequences
%     Nn = GA.Population;
%     ids = 1:Nn;
    % Get only top %
	ids = GA.GetTopPop(GA.Fittest(1));
    Nn = GA.Fittest(1);

    if strcmp(type, 'seq')
%         data = GA.Seqs(:,:,Gen);
        data = GA.Seqs(ids,:,Gen);
        
        % Normalize data
        data = bsxfun(@rdivide, ...
            bsxfun(@minus, data, GA.Gen.Range(1,:)), ...
            (GA.Gen.Range(2,:)-GA.Gen.Range(1,:)));
    else
        data = GA.Fit(ids,:,Gen);
    end
    
    dist = zeros(Nn,1);
    for i = 1:Nn
        % Find closest neighbor
        dist_mat = data - repmat(data(i,:),Nn,1);
        dist_vec = sqrt(diag(dist_mat*dist_mat'));
        
        % Remove self
        dist_vec(i) = [];
        dist(i) = min(dist_vec);
    end
end

