function [Lines] = ConnectVectors( Lines )
    N=size(Lines,2);

    for d=1:N-1
        L=length(find(Lines(:,d)~=0));
        if L==10
            % To shorten the computation time of calculating
            % all the permutations of a 10 elements vector,
            % we remove the 2 eigenvalues that are closest
            % to zero (there's a zero eigenvalue in the PM)

            CurVec=[Lines(1:L,d), (1:L)']; % Vector that holds the current points
            NextVec=[Lines(1:L,d+1), (1:L)']; % Vector that holds the next set of points

            CurVecS=sortrows(CurVec,1);
            NextVecS=sortrows(NextVec,1);

            % Grab the last 8
            CurVec=CurVecS(3:end,1);
            NextVec=NextVecS(3:end,1);

            L=8;
            if length(CurVec)~=length(NextVec)
                A=1;
            end
        else
            CurVec=Lines(1:L,d); % Vector that holds the current points
            NextVec=Lines(1:L,d+1); % Vector that holds the next set of points
        end

        % Calculate all possible permutations of the next points
        % This takes a LONG LONG LONG while when M is large
        PermVec=perms(NextVec);

        % Calculate distance from current vector to all permutations
        Dist=abs(PermVec-repmat(CurVec',size(PermVec,1),1));

        % Find the most distant point for each permutation
        MaxDist=max(Dist,[],2);

        % Find the permutation where the max. distance is the smallest
        MinPerm=find(MaxDist==min(MaxDist));

        MP=length(MinPerm);
        if MP==1
            NextVec=PermVec(MinPerm,:)';
            if L~=8
                Lines(1:L,d+1)=NextVec;
            else
                Lines(CurVecS(1:2,2),d+1)=Lines(NextVecS(1:2,2),d+1);
                Lines(CurVecS(3:end,2),d+1)=NextVec;
            end
        else            
            % Too many best permutations found
            % We'll select the one with the minimal overall distance
            BestPerms=PermVec(MinPerm,:);

            % Calculate absoulte distance from CurVec to each candidate
            Dist=abs(BestPerms-repmat(CurVec',MP,1));
            AbsDist=diag(Dist*Dist');

            MinPerm=find(AbsDist==min(AbsDist));
            if length(MinPerm)>1
                % screw it, choose the first one
                NextVec=BestPerms(MinPerm(1),:)';
            else
                NextVec=BestPerms(MinPerm,:)';
            end

            if L~=8
                Lines(1:L,d+1)=NextVec;
            else
                Lines(CurVecS(1:2,2),d+1)=Lines(NextVecS(1:2,2),d+1);
                Lines(CurVecS(3:end,2),d+1)=NextVec;
            end
        end
    end
end
