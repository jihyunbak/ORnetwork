function ndat_cell = plotMultiLevelTree(gall,nodenames,mycolors,colorByLevel,opts)
% plotMultiLevelTree: Plots a multi-level tree with hard clustering.
% 
% INPUT:
%   - gall: [N L] array, N: # nodes, L: # grouping levels.
%           each [N 1] column is a list of group indices at the level.
%   - nodenames: [N 1] cell array of strings.
%   - mycolors: a [K 3] array, with a list of RGB colors.
%               K should be no less than the largest group index
%               at the target level (specified by colorByLevel).
%   - colorByLevel: specifies to which level of grouping color should be
%                   applied. corresponds to a column index for gall.
%   - opts: various options for plot. See content for details.
% OUTPUT:
%   - ndat_cell: a cell array of nodes, with group data.

% Copyright 2018 Ji Hyun Bak
% ------------------------------------------------------------------------


%% unpack input data

numNodes = size(gall,1);
numLevels = size(gall,2);

% unpack options (if not specified, set to default value)
fsz = getFromStruct(opts,'fontsize',14);
msz = getFromStruct(opts,'markersize',10);
lwd = getFromStruct(opts,'linewidth',1);
lcol = getFromStruct(opts,'linecolor',0.5*[1 1 1]);
colorChoiceRule = getFromStruct(opts,'colorchoicerule','first'); % 'first' or 'last'
plotAxis = getFromStruct(opts,'plotaxis',1); % 1 (x) or 2 (y)
xreverse = getFromStruct(opts,'xreverse',false);
yreverse = getFromStruct(opts,'yreverse',false);
textgap = getFromStruct(opts,'textgap',0.1); 

%% calculate node position

ndat_cell = cell(numLevels,1); % node data
for nl = 1:numLevels
    myglist = gall(:,nl);
    myguniq = unique(myglist,'stable');
    if(nl==1)
        myx = (1:numel(myguniq))';
    else
        mygprev = ndat_cell{nl-1}.g;
        myxprev = ndat_cell{nl-1}.x;
        
        myx = zeros(numel(myguniq),1);
        for g = myguniq(:)'
            gprev = unique(gall(gall(:,nl)==g,nl-1));
            myx(myguniq==g) = mean(myxprev(ismember(mygprev,gprev)));
        end
    end
    mygnext = zeros(numel(myguniq),1);
    mygcol = zeros(numel(myguniq),1);
    for i = 1:numel(myguniq)
        myg = myguniq(i);
        if(nl<numLevels)
            mygnextuniq = unique(gall(gall(:,nl)==myg,nl+1)); % normally should be only 1
            if(numel(mygnextuniq)>1)
                error('non-hierarchical grouping');
            end
            mygnext(i) = mygnextuniq;
        end
        mygcol_all = find(gall(:,nl)==myg,1,colorChoiceRule); % may not be unique! follow colorChoiceRule
        mygcol(i) = gall(mygcol_all,colorByLevel); % node color index
    end
    
    % save node data
    ndat_cell{nl} = struct('g',myguniq,'x',myx,'gnext',mygnext,'gcol',mygcol);
    
end

%% plot multi-level tree graph

% figure(21)
% clf;
for nl = 1:numLevels

    myndat = ndat_cell{nl};
    myguniq = myndat.g;
    myx = myndat.x;
    mygnext = myndat.gnext;
    mygcol = myndat.gcol;
    
    for i = 1:numel(myguniq)
        % edges
        if(nl<numLevels)
            myndat_next = ndat_cell{nl+1};
            myxnext = myndat_next.x(myndat_next.g==mygnext(i));
            if(plotAxis==1)
                myplotx = [myx(i) myxnext];
                myploty = [nl nl+1];
            elseif(plotAxis==2)
                myplotx = [nl nl+1];
                myploty = [myx(i) myxnext];
            else
                error('unknown plot axis option');
            end
                plot(myplotx,myploty,'-','color',lcol,'linewidth',lwd)
            hold on
        end
        % nodes
        if(plotAxis==1)
            plot(myx(i),nl,'ko','markerfacecolor',mycolors(mygcol(i),:),'markersize',msz)
            if(nl==1 && ~isempty(nodenames{i}))
                if(yreverse)
                    text(myx(i),nl-textgap,nodenames{i},'verticalalignment','middle','rotation',90,'fontsize',fsz)
                else
                    text(myx(i),nl-textgap,nodenames{i},'verticalalignment','middle','rotation',-90,'fontsize',fsz)
                end
            end
        elseif(plotAxis==2)
            plot(nl,myx(i),'ko','markerfacecolor',mycolors(mygcol(i),:),'markersize',msz)
            if(nl==1 && ~isempty(nodenames{i}))
                if(xreverse)
                    text(nl-textgap,myx(i),nodenames{i},'horizontalAlignment','left','fontsize',fsz)
                else
                    text(nl-textgap,myx(i),nodenames{i},'horizontalAlignment','right','fontsize',fsz)
                end
            end
        else
            error('unknown plot axis option');
        end
        hold on
    end
end
hold off
if(yreverse)
    set(gca,'ydir','reverse')
end
if(xreverse)
    set(gca,'xdir','reverse')
end
axis tight
% xlim([0.5 numNodes+0.5])
% ylim([0.5 numLevels+0.5])


end

function myval = getFromStruct(mystruct,myfield,defaultval)
% inherit specified field from a struct if the field exists,
% otherwise set to default value
%
% INPUTS:
% - mystruct: a struct variable
% - myfield: a string (supposed to be the field name)
% - defaultval: default value to set when field is nonexisting

if(isfield(mystruct,myfield))
    myval = mystruct.(myfield);
else
    myval = defaultval; % set default value
end

end
