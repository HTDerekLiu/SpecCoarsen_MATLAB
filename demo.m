clc; clear all; close all;
addpath('./utils/')
% addpath('path/to/gptoolbox') (gptoolbox can be found in "https://github.com/alecjacobson/gptoolbox")

% parameters
numEig = 100; % number of eigenfunctions to preserve
numNc = 500; % number of coarse points
lr = 2e-2; % learning rate
decayIter = 1; % learning rate decay iterations (it is optional, just for fine tune the result)
stallIter = 5; % stalling iteration 

% read mesh
[V,F] = readOBJ('./bunny.obj');
V = normalizeUnitArea(V,F);

% construct an initial operator and mass matrix
L = -cotmatrix(V,F);
M = massmatrix(V,F);

% algebraic coarsening
% note: this matlab implementation does not implement the sparse gradient in the 
% appendix A of "Spectral Coarsening of Geometric Operators" [Liu et al. 2019].
% Thus it would be much slower than C++ implementation.
[Lc, Mc, G, P, Cpt] = algebraicCoarsening(L, M, numNc, numEig, ...
    'lr', lr, 'decayIter', decayIter, 'stallIter', stallIter);
 
% visualize functional map
[~, eVecc] = eigsReal(Lc, Mc, numEig);
[~, eVec] = eigsReal(L, M, numEig);
fMap = eVecc' * Mc * P * eVec;
figure(1)
plotFMap(fMap)
title('functional map image')

% visualize one eigenfunctions (this may have sign flip)
figure(2)
subplot(1,2,1)
plotMesh(V,F,eVec(:,10))
subplot(1,2,2)
scatter3(V(Cpt,1),V(Cpt,2),V(Cpt,3), 30, eVecc(:,10), 'filled')
axis equal off
title('visualize one eigen functioon')

% visualize root nodes
figure(3)
plotMesh(V,F)
hold on
scatter3(V(Cpt,1),V(Cpt,2),V(Cpt,3),20,'filled')
title('visualize root nodes')
