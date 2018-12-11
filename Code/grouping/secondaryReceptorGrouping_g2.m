function [rgroup2,G12map] = secondaryReceptorGrouping_g2(rgroup1,chi_g1g1,ovlpcut)
% secondaryReceptorGrouping_g2:
% Performs a greedy pairwise merging of the primary receptors groups
% at a specified threshold.
% 
% INPUT:
%   - rgroup1: [N 1] vector of primary group indices, N: # receptors
%   - chi_g1g1: inter-group overlaps between primary receptor groups 
%   - ovlpcut: cutoff for pairwise merging
% OUTPUT:
%   - rgroup2: [N 1] vector of secondayr group indices, ranked by size
%   - G12map: the merging operation that maps rgroup1 to rgroup2

% Copyright 2018 Ji Hyun Bak
% ------------------------------------------------------------------------

% receptor grouping g2: merge g1 groups by hard overlap cutoff
G12map_raw = pairwiseMerge_hardCutoff(chi_g1g1,ovlpcut);
rgroup2raw = G12map_raw(rgroup1);

% re-number group indices
rgroup2 = renumberGroupsBySize(rgroup2raw); % sort by group size (total # elements)

% summarize the merging operation (after re-numbering)
rg12_pairs = unique([rgroup1(:) rgroup2(:)],'rows'); 
G12map = zeros(max(rgroup1),1); % transition operator
for myg = 1:max(rgroup1)
    G12map(myg) = rg12_pairs(rg12_pairs(:,1)==myg,2);
end

end

% ------------------------------------------------------------------------

function mygroup = pairwiseMerge_hardCutoff(pairmat,mycutoff)
% greedy pairwise merging with pairwise scores and a hard cutoff

% apply hard cutoff
[j1,j2] = find(pairmat>=mycutoff);
myedges = [j1 j2]; % all pairs to be merged
myedges = myedges(j1<j2,:); % remove duplicates
mynodes = (1:size(pairmat,1))';

% greedy merging
mygroup = mynodes; % the "group leader" index (earliest node)
for j = size(myedges,1):-1:1
    % group the connected nodes (and their group members) together
    mygroup(ismember(mygroup,mygroup(myedges(j,:)))) = min(mygroup(myedges(j,:)));
end

% check merge
grpcheck = bsxfun(@eq,mygroup(myedges(:,1)),mygroup(myedges(:,2)));
if(~all(grpcheck))
    error('pairwiseMerge: grouping error');
end

end
