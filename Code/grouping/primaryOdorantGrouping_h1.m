function lgroup1 = primaryOdorantGrouping_h1(myBIM,rgroup1)
% primaryOdorantGrouping_h1:
% Assigns each odorant to the most-interacting primary receptor group, 
% based on the receptor code overlap.
% 
% INPUT
%   - myBIM: [N M] binary interaction matrix
%     each column of myBIM is a receptor code (length N = # receptors)
%     number of columns M = # odorants
%   - rgroup1: [N 1] vector of recepetor group indices, integer values
% OUTPUT
%   - lgroup1: [M 1] vector of odorant group indices,
%              values inherited from rgroup1

% Copyright 2018 Ji Hyun Bak
% ------------------------------------------------------------------------



% unpack input: collect interacting pairs
[pairR,pairL] = find(myBIM>0);
numL = size(myBIM,2);

% set up receptor bases (inherit from receptor groups)
allgroups = cell(max(rgroup1),1);
for rg = 1:max(rgroup1)
    myRs = find(rgroup1==rg);
    allgroups{rg} = myRs;
end

% assign each odorant to the most-overlapping group
lgroup1 = zeros(numL,1);
for jL = 1:numel(lgroup1)
    myRs = pairR(pairL==jL); % receptor code for this odorant
    ovlpvec = calculateAllOverlaps(myRs,allgroups); % calculate all overlaps
    maxovlp = max(ovlpvec);
    newg = find(ovlpvec==maxovlp,1,'first'); % pick highest-rank group
    lgroup1(jL) = newg;
end

end

% ------------------------------------------------------------------------

function ovlpvec = calculateAllOverlaps(myRs,allgroups)

% set up normalized overlap function
getoverlap = @(yj,ybase) numel(intersect(ybase,yj))/numel(yj); % unweighted

% calculate all overlaps
ovlpvec = zeros(numel(allgroups),1); % dimension-reduced vector
for ngrp = 1:numel(allgroups)
    thegroup = allgroups{ngrp};
    theovlp = getoverlap(myRs,thegroup); % unweighted overlap
    ovlpvec(ngrp) = theovlp;
end

end
