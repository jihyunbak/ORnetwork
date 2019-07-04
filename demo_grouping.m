% demo_grouping.m
% 
% Demo script for the grouping algorithms, 
% applied to the binary interactions between odorants and receptors.

% Copyright 2018 Ji Hyun Bak

%% initialize

clear all;
clc;

addpath(genpath('Code')); % add path to custom code


%% load pairwise interaction data

% input directory
mydir = 'Data/database/';

% all nodes
nodefilename = [mydir,'ORnet_node_table.csv'];
nodetab = readtable(nodefilename);

% all edges
edgefilename = [mydir,'ORnet_edge_table.csv'];
edgetab = readtable(edgefilename);


%% the interaction matrix


% ===== construct the interaction matrix

% unique indices for receptors and odorants, for all interacting pairs
iR0 = edgetab.Ridx; % receptors (R)
iL0 = edgetab.Lidx; % odorants (L; ligands)

% interaction matrix
pairidxmat = sparse(iR0,iL0,(1:numel(iR0))'); % with edge index
myBIM = (pairidxmat>0); % binary interaction matrix

[numR,numL] = size(myBIM); % [N M], or [#receptors #odorants]


% ===== plot 

% open a figure window
clf;

% force portrait view
fh = gcf;
mypos = fh.Position;
if(mypos(4)<mypos(3)) % if landscape
    mypos(4) = mypos(3)*1.2; % make portrait type
    mypos(2) = mypos(2)-mypos(3)*0.2; % avoid overflow
end
set(gcf,'Position',mypos)

% plot interaction matrix
subplot(2,2,1)
imagesc(myBIM')
xlabel('receptors')
ylabel('odorants')
title('unsorted interaction matrix')
colormap(gca,[1 1 1; 0 0 0])
set(gca,'ydir','normal')

drawnow;


%% the co-activation matrix of receptors

% ===== construct the co-activation matrix of receptors

% co-activation matrix of receptors
rcomat0 = double(double(myBIM)*double(myBIM)'>0); % [N N] binary matrix


% ===== plot

subplot(2,2,3)
imagesc(rcomat0)
axis square
xlabel('receptors')
ylabel('receptors')
title('unsorted co-activation matrix')
colormap(gca,[1 1 1; 0 0 0])
set(gca,'ydir','normal')

drawnow;


%% ranking odorants and receptors separately (g0 and h0)

% ===== ranking

% odorant rank h0: by degree
degL = full(sum(myBIM,1)'); % # receptors activated by each odorant
h0rank = rankBy(degL,'descend'); % 12/11/2018
[~,ilsrt0] = sort(h0rank);

% receptor rank g0: by degree of co-activation
degRcoact = full(sum(rcomat0,2)); % degree in (binarized) coactivation graph
g0rank = rankBy(degRcoact,'descend'); % 12/11/2018
[~,irsrt0] = sort(g0rank);


% ===== plot sorted results

% plot interaction matrix
subplot(2,2,2)
imagesc(myBIM(irsrt0,ilsrt0)')
xlabel('receptors')
ylabel('odorants')
title('sorted by g0/h0')
colormap(gca,[1 1 1; 0 0 0])
set(gca,'ydir','normal')

% plot co-activaion matrix
subplot(2,2,4)
imagesc(rcomat0(irsrt0,irsrt0))
axis square
xlabel('receptors')
ylabel('receptors')
title('sorted by g0')
colormap(gca,[1 1 1; 0 0 0])
set(gca,'ydir','normal')

drawnow;


%% primary grouping


% ===== primary grouping

% primary receptor grouping g1
rgroup1 = primaryReceptorGrouping_g1(myBIM);
[~,irsrt1] = sortrows([g0rank rgroup1],[2 1]); % sort by g1 then g0

% primary odorant grouping h1
lgroup1 = primaryOdorantGrouping_h1(myBIM,rgroup1);
[~,ilsrt1] = sortrows([h0rank lgroup1],[2 1]); % sort by h1 then h0



% ===== sort by primary grouping, and plot

% sorted interaction matrix
subplot(2,2,1)
imagesc(myBIM(irsrt1,ilsrt1)')
xlabel('receptors')
ylabel('odorants')
title('sorted by g1/h1')
colormap(gca,[1 1 1;0 0 0])
set(gca,'ydir','normal')

% sorted co-activation matrix
subplot(2,2,3)
imagesc(rcomat0(irsrt1,irsrt1))
axis square
xlabel('receptors')
ylabel('receptors')
title('sorted by g1')
colormap(gca,[1 1 1; zeros(max(rgroup1),3)])
set(gca,'ydir','normal')

drawnow;


%% secondary grouping


% ===== secondary grouping

% get inter-group overlap between g1 receptor groups
chi_g1g1 = getInterGroupOverlap(rcomat0,rgroup1);

% secondary receptor grouping g2 
ovlpcut = 0.15; % (at fixed optimal threshold)
[rgroup2,G12map] = secondaryReceptorGrouping_g2(rgroup1,chi_g1g1,ovlpcut); 
[~,irsrt2] = sortrows([g0rank rgroup1 rgroup2],[3 2 1]); % sort by g2 then g1 then g0
    
% secondary odorant grouping h2
lgroup2 = secondaryOdorantGrouping_h2(lgroup1,G12map);
[~,ilsrt2] = sortrows([h0rank lgroup1 lgroup2],[3 2 1]); % sort by h2 then h1 then h0



% ===== set up colors for the secondary groups

% colors for the 6 largest groups
selectColors = [... 
    1 0 0; ... % red
    1 0.5 0; ... % orange
    0.3 0.9 0.2; ... % dark green
    0.2 0.8 1; ... % cyan
    0 0 1; ... % blue
    2/3 0 1 ... % purple
    ];

% coloring by g2 index: 
% the first (largest) 6 groups only; all other groups in gray
mycmap_by_g2 = [selectColors; 0.7*ones((max(rgroup2)-6),3)]; % (trimmed at group 6)



% ===== plot

% interaction matrix, sorted and labeled by secondary receptor groups
myBIM_srt2 = bsxfun(@times,rgroup2(irsrt2),myBIM(irsrt2,ilsrt2));

subplot(2,2,2)
imagesc(myBIM_srt2')
xlabel('receptors')
ylabel('odorants')
title('sorted by g2/h2')
colormap(gca,[1 1 1; mycmap_by_g2])
set(gca,'ydir','normal')


% co-activation matrix, sorted and labeled by secondary groups
rcomat_by_g2 = bsxfun(@times,...
    bsxfun(@max,rgroup2(irsrt2),rgroup2(irsrt2)'),... % color with smaller group when overlapping
    rcomat0(irsrt2,irsrt2));

subplot(2,2,4)
imagesc(rcomat_by_g2)
axis square
xlabel('receptors')
ylabel('receptors')
title('sorted by g2')
colormap(gca,[1 1 1; mycmap_by_g2])
set(gca,'ydir','normal')


%% plotting all together: reproducing Fig 4 in the manuscript


% ===== prepare for odorant plotting

% odorant names (sort by secondary grouping)
odorant_tab = nodetab(nodetab.NodeType==2,:);
odorant_names_srt_noquote = strrep(odorant_tab.LongName(ilsrt2),'"',''); % sorted to match plot

% receptor names (sort by secondary grouping)
receptor_tab = nodetab(nodetab.NodeType==1,:);
receptor_names_srt = receptor_tab.LongName(irsrt2); % sorted to match plot

% choose if you want odorant/receptor names on plot
printOdorantNames = true;
printReceptorNames = false;
if(printOdorantNames)
    Lnodenames = odorant_names_srt_noquote;
else
    Lnodenames = cell(numL,1);
end
if(printReceptorNames)
    Rnodenames = receptor_names_srt;
else
    Rnodenames = cell(numR,1);
end


% ===== plot, using custom function plotMultiLevelTree

w0 = 0.05;
h0 = 0.05;
wgap = 0.05;
hgap = 0.05;
% for the activation matrix
wA = 0.6;
hA = 0.75;
% for the trees
wT = 1-wA-2*w0-wgap;
hT = 1-hA-2*h0-hgap;

% options for multi-level tree plots
fsz = 10;
target_lev = 3; % target for group coloring (secondary group)

clf;

% --- activation matrix
axes('position',[w0+wT+wgap h0+hT+hgap wA hA])
imagesc(myBIM_srt2')
colormap(gca,[1 1 1;mycmap_by_g2]) % white for 0
set(gca,'ydir','normal','linewidth',1,'box','off')
xlabel('receptors (sorted)')
ylabel('odorants (sorted)')
xlim([0 numR])
ylim([0 numL])

% --- odorant tree
axes('position',[w0 h0+hT+hgap wT hA])
opts = struct('plotaxis',2,'yreverse',false,'xreverse',true,'fontsize',fsz,'textgap',0.3);
plotMultiLevelTree([(1:numL)' lgroup1(ilsrt2) lgroup2(ilsrt2)],...
    Lnodenames,mycmap_by_g2,target_lev,opts);
xlim([-5 2])
axis off

% --- receptor tree
axes('position',[w0+wT+wgap h0 wA hT])
opts = struct('yreverse',true,'fontsize',fsz);
plotMultiLevelTree([(1:numR)' rgroup1(irsrt2) rgroup2(irsrt2)],...
    Rnodenames,mycmap_by_g2,target_lev,opts);
ylim([1 4])
axis off

