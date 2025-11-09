% DTucker: Given a dense tensor, decompose the tensor into factor matrices and core tensor
% 
% Parameter
% X: input tensor which has dense form
% J: target rank list
% tol: tolerance on difference (stop criterion)
% maxiters: maximum number of iterations 
% init: wheter to use initialization phase (=1) or not (=0)
%
% Return
% core: the core tensor of tucker decomposition
% factor matrices: factor matrices of tucker decomposition
%

function [core, factors] = DTucker(X, J, tol, maxiters, init)

% decomp, init, update, num_iter
% experim = [0,0,0,0];

fitchangetol = tol;

N = ndims(X);
normX = norm(X);

X = tensor(X);


% tstart=tic;
% tensorSize = size(X);
[~,idx] = sort(size(X), 'descend');

% 定义要调换的元素索引
i = 3;
j = 4;
% 调换元素顺序
temp = idx(i);
idx(i) = idx(j);
idx(j) = temp;


X = permute(X, idx);
Jtmp = J;
for rank_permute = 1:N
    J(rank_permute) = Jtmp(idx(rank_permute));
end

dims = size(X);
Jsort = sort(J, 'descend');
[decompU, decompS, decompV] = Ddecomp(X, Jsort(2));

% tElapsed = toc(tstart);
% experim(1) = tElapsed;

Ainit = cell(1, N);
for i=1:N
    Ainit{i} = randn(dims(i), J(i));
end

% tstart=tic;
% init;
if init == 1
    Ainit = Dinit(decompU, decompS, decompV, dims, J);
end
% tElapsed = toc(tstart);
% experim(2) = tElapsed;

% test = Ainit{3}.'*Ainit{3};

fit = 0;
A = Ainit;
num_iter = 0;
% tstart = tic;
for iter=1:maxiters
%     tic;
    fitold = fit;
    [A, G] = Dupdate(decompU, decompS, decompV, A, dims); 
    % compute error
    normresidual = sqrt( normX^2 - norm(tensor(G))^2 );
    fit = 1 - (normresidual / normX); %fraction explained by model
    fitchange = abs(fitold - fit);

%     fprintf(' Iter %2d: fit = %e fitdelta = %7.1e\n', iter, fit, fitchange);
%     toc;
    % Check for convergence
    if (iter > 1) && (fitchange < fitchangetol)
        num_iter = iter;
        break;
    end

end

if num_iter == 0
    num_iter = maxiters;
end

% tElapsed = toc(tstart);
% experim(3) = tElapsed;
% experim(4) = num_iter;

% fprintf(' Running time of approximation phase is  %4f seconds\n', experim(1));
% fprintf(' Running time of initialization phase is  %4f seconds\n', experim(2));
% fprintf(' Running time of iteration phase is  %4f seconds\n', experim(3));
% fprintf(' Number of iteration is  %2d\n', experim(4));

% X = ipermute(X, idx);
G = ipermute(G, idx);
Atmp = A;
for a=1:N
   A{idx(a)} = Atmp{a}; 
end
clear Atmp;

factors = A;
core = G;
end
