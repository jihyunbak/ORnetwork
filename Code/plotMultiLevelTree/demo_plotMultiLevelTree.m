% demo_plotMultiLevelTree.m

% demo script for custom function plotMultiLevelTree.m
% 12/11/2018 JHB

%% initialize

clear all;


%% generate example

nodenames = {'A11';'A12';'A13';'A21';'A22';'B11';'B12';'B21';'B22'};

g1 = 1:9; % original nodes
g2 = [1 1 1 2 2 3 3 4 4];
g3 = [1 1 1 1 1 2 2 2 2];
gall = [g1(:) g2(:) g3(:)];

mycolors = colormap('prism');
colorByLevel = 2; % color nodes by group index at this level


%% test plot

% set various options (see inside functions for details)
opts = struct('fontsize',12,'markersize',15,'linewidth',2,'plotaxis',2,'yreverse',true);

% plot
clf;
plotMultiLevelTree(gall,nodenames,mycolors,colorByLevel,opts)
axis off
