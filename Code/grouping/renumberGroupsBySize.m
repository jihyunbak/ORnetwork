function grouping = renumberGroupsBySize(grouping_old)
% renumberGroupsBySize: Re-numbers group indices from the largest group,
% while keeping the original group index order in case of ties.
% 
% INPUT: grouping_old: a vector of grouping indices
% OUTPUT: grouping: a vector of matching size, with re-numbered indices

% Copyright 2018 Ji Hyun Bak
% ------------------------------------------------------------------------

grouping = -grouping_old; % negative sign to indicate aux idx

% detect group sizes
[grpIdxNeg,~,ib] = unique(grouping,'rows');
nodeCnt = accumarray(ib,1);

% sort by group sizes, in descending order
mysrt = sortrows([grpIdxNeg nodeCnt],[-2 -1]); % also keep original grp Idx order (-1)

% update numbering
for ng = 1:size(mysrt,1)
    grouping(grouping==mysrt(ng,1)) = ng; % largest group first
end

end
