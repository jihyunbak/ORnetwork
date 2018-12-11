function myrank = rankBy(V,varargin)
% rankBy: Rank a list of elements V by their values, 
% in the order specified in sortOrder.
% 
% INPUT
%   - V: a list of values, [N 1] vector
%   - sortOrder [optional]: order of sorting, 
%                           either 'ascend' (default) or 'descend'.
% OUTPUT
%   - myrank: a list of integers (the ranks), spanning the range 1:N.

% Copyright 2018 Ji Hyun Bak
% ------------------------------------------------------------------------

sortOrder = 'ascend'; % default
if(nargin>1)
    sortOrder = varargin{1};
end

idx0 = (1:numel(V))'; % original order

[~,idx_aux1] = sort(V,sortOrder); % 12/11/2018 test
[~,idx_aux2,~] = unique(idx_aux1); % inverse sort
myrank = idx0(idx_aux2); % h0 rank

end
