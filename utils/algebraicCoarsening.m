function [Lc, Mc, G, P, Cpt] = algebraicCoarsening(L, M, m, numEig, varargin)
% ALGEBRAICCOARSENING simplifies a given PSD operator L with size n-by-n 
% to m-by-m using the method of 
%   Spectral Coarsening of Geometric Operators [Liu et al. 2019]
%
% [Lc, Mc, G, P, Cpt] = algebraicCoarsening(L, M, m, numEig)
%
% Inputs:
%   L n x n PSD matrix of a differential operator 
%    (diagonal positive, off-diagonal mostly negative)
%   M n x n diagonal matrix of variable masses (e.g. vertex areas)
%   m a integer number controls the size of coarsened matrix
%   numEig a integer number controls how many eigenvectors to preserve
% Outputs:
%   Lc a m x m matrix of a coarsened differential operator
%   Mc a m x m matrix of a coarsened mass matrix
%   G the interpolation operator where Lc = G' * L * G;
%   P the projection operator from dense to the coarsened domain Uc = P * U
%   Cpt the indices of root nodes 

% default parameter values
lr = 2e-2; % learning rate (step size) for the gradient descent
lrDecayIter = 2; % learning rate decay iterations (optional)
maxIter = 1000; % the maximum iteration (usually won't reach)
stallIter = 5; % stall iteration

% process user input parameters
params_to_variables = containers.Map(...
    {'lr', 'decayIter', 'stallIter', 'maxIter'}, ...
    {'lr', 'lrDecayIter', 'stallIter', 'maxIter'});
v = 1;
while v <= numel(varargin)
    param_name = varargin{v};
    if isKey(params_to_variables,param_name)
        assert(v+1<=numel(varargin));
        v = v+1;
        % Trick: use feval on anonymous function to use assignin to this workspace 
        feval(@()assignin('caller',params_to_variables(param_name),varargin{v}));
    else
        error('Unsupported parameter: %s',varargin{v});
    end
    v=v+1;
end

%% Combinatorial Coarsening
invM = diag(diag(M).^-1);
SoC = strengthOfConnections(L, M);

% clustering
seedsIdx = randperm(size(L,1));
seedsIdx = seedsIdx(1:m);
clusterAssignment = graphKMediods(SoC, seedsIdx);

% construct P (projection) and K (assignment)
[Cpt, ~, idxToCpt] = unique(clusterAssignment);
P = sparse([1:length(Cpt)], Cpt, ones(length(Cpt),1), ...
    length(Cpt), size(L,1), length(Cpt));
K = sparse(idxToCpt, [1:size(L,1)], ones(length(idxToCpt),1), ...
    length(Cpt), size(L,1), length(idxToCpt));

%% Operator Optimization
Mc = K * M * K'; % coarse mass matrix
invMc = diag(diag(Mc).^-1);

SL = double(L~=0);
Ac = double(K * SL * K' ~= 0);
SG = double(K' * Ac ~= 0);
% SLc = double(SG' * SL * SG ~= 0);

% eigen vector as test vectors
[~, U] = eigsReal(L, M, numEig);

G = K'; % initial G

% energy function/gradient precomputation
A = P * invM * L * U;
B = P * U;
BB = B*B';
AB = A*B';

% projection to null space precomputation
U0 = U(:,1);
ZIdx = find(SG);
[rIdx, cIdx] = find(SG);
Z = sparse(ZIdx, [1:length(ZIdx)], ones(length(ZIdx),1), size(G,1)*size(G,2), length(ZIdx));
g = G(ZIdx);
Aproj = kron((P*U0)',speye(size(G,1))) * Z;
AAA = Aproj'*inv(Aproj*Aproj');

% stopping criteria 
objValOld = 1e10;
stop = 0;
lrDecayCount = 0;
decayRatio = 0.5;

% NADAM
timeStep = 0;
beta1 = 0.9;
beta2 = 0.9;
eps = 1e-8;
mt = sparse(size(G,1), size(G,2));
nt = sparse(size(G,1), size(G,2));

% function for computing objective function
objFunc = @(params) objFuncTemplate(params,A, B, L, Mc, invMc);
errorHis = [];

for iter = 1:maxIter
    % precomputation
    objVal = objFunc(G);
    errorHis = [errorHis objVal];

    % compute gradient 
    % Note that MATLAB implementation DOES NOT exploit the sparse gradient 
    % mentioned in the appendix A of [Liu et al. 2019] because looping 
    % over the indices is slow in MATLAB. Therefore, this is much slower 
    % than c++ implementation
    LG = L*G;
    BBGLGinvMc = BB*(G'*(LG*invMc));
    grad = LG * (- AB' - AB + BBGLGinvMc + BBGLGinvMc') ;
    grad = grad .* SG;

    % NADAM
    timeStep = timeStep + 1;
    mt = beta1*mt + (1-beta1)*grad;
    nt = beta2*nt + (1-beta2)*(grad.^2);
    mt_hat = mt / (1-beta1^timeStep);
    nt_hat = nt / (1-beta2^timeStep);
    Ngrad = 1./(sqrt(nt_hat) + eps) .* (beta1*mt_hat + (1-beta1)*grad/(1-beta1^timeStep));

    % update G
    G = G - lr * Ngrad;
    g1 = G(ZIdx);
    g = g1 - AAA * (Aproj*g1-U0(:));
    G = sparse(rIdx, cIdx, g, size(SG,1), size(SG,2));
    
    if (max(abs(G * P * U0 - U0)) > 1e-6)
        fprintf('projection failed\n')
    end

    % print progress
    if mod(iter,10) == 0
        fprintf('iter %i, cost %f\n', iter, objVal);
    end

    % stopping criteria (if the energy stop decreasing for a few iterations)
    if objValOld < objVal
        stop = stop + 1;
        if stop > stallIter
            if lrDecayCount < lrDecayIter % not decrease lr yet
                lr = lr * decayRatio;
                stop = 0;
                lrDecayCount = lrDecayCount + 1;
                fprintf('derease learning rate to: %f\n', lr);
            else
                fprintf('iter %i, cost %f\n', iter, objValOld);
                break;
            end
        end
    else
        stop = 0;
        objValOld = objVal;
        Gbest = G;
    end
end
G = Gbest;
Lc = G' * L * G;
Lc = (Lc + Lc')/2; % sometimes a slightly non-symmetric matrix fails 'eigs'