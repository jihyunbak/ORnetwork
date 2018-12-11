function rgroup1 = primaryReceptorGrouping_g1(BIM_input)
% primaryReceptorGrouping_g1:
% Makes non-overlapping partitions of co-activated receptor *cliques*
% associated to shared hub odorants.
% 
% INPUT:  BIM_input, [N M] binary interaction matrix
%                    N: number of receptors, M: number of odorants
% OUTPUT: rgroup, [N 1] vector of group indices, with integer values

% Copyright 2018 Ji Hyun Bak
% ------------------------------------------------------------------------

%% unpack input 

myBIM = (BIM_input>0); % force binary
numR = size(myBIM,1); % number of receptors (rows)


%% primary receptor grouping

% sort the interaction matrix by odorant degree
degL = full(sum(myBIM,1)); % odorant degree: # receptors interacting
[~,ilsrt0] = sort(degL,'descend'); % this corresponds to the h0 rank
myBIM_srtdeg = full(myBIM(:,ilsrt0));

% group by highest-degree odorant partners
rgroup1raw = zeros(numR,1);
for nrs = 1:numR
    rgroup1raw(nrs) = find(myBIM_srtdeg(nrs,:)>0,1,'first');
end

% re-number group indices
rgroup1 = renumberGroupsBySize(rgroup1raw); % sort by group size

end
