function chi_g1g1 = getInterGroupOverlap(pairmat,grpvec)
% getInterGroupOverlap:
% Calculates the inter-group overlap in a matrix of pairwise connections
% given a grouping that maps each element to a group index.
% 
% INPUT:
%   - pairmat: [N N] symmetric matrix with pairwise connetions.
%   - grpvec: [N 1] vector of group indices, with integer values.
% OUTPUT:
%   - chi_g1g1: [K K] symmetric matrix of inter-group overlap,
%                where K is the largest group index in grpvec.

% Copyright 2018 Ji Hyun Bak
% ------------------------------------------------------------------------

% function for inter-group overlap calculation
getChi = @(ig1,ig2) getChi_offblock(pairmat,ig1,ig2,grpvec);

% inter-group overlaps
chi_g1g1 = NaN(max(grpvec),max(grpvec));
for ig1 = 1:max(grpvec)
    for ig2 = ig1:max(grpvec)
        mychi = getChi(ig1,ig2);
        chi_g1g1(ig1,ig2) = mychi;
        chi_g1g1(ig2,ig1) = mychi;
    end
end

end

% ------------------------------------------------------------------------

function mychi = getChi_offblock(pairmat,ig1,ig2,grpvec)
% count the fraction of non-zero elements in the matrix (binary limit)
% across different groups

% across-group (off-diagonal-block) submatrix
mybox = pairmat(grpvec==ig1,grpvec==ig2); 

% overlap fraction
mychi = sum(mybox(:)>0)/numel(mybox); 

end
