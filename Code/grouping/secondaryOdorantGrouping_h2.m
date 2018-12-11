function lgroup2 = secondaryOdorantGrouping_h2(lgroup1,G12map)
% secondaryOdorantGrouping_h2:
% Seconday grouping of odorants, which is a simple propagation of the
% merging operation determined from the receptor groups.
% > lgroup2 = G12*lgroup1. 

% Copyright 2018 Ji Hyun Bak
% ------------------------------------------------------------------------

% -- simply apply the merging operation
lgroup2 = G12map(lgroup1);

