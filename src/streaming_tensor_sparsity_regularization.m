
function [As, Ps_update, core, Pcore_update, Qcore_update, alpha ]  =  ITSReg_update(D, As, Ps, core, Pcore, Qcore,beta,Par,maxIter)

%% 
N = length(As)+1;
Nlist = 1:N-1;
sizeD          = size(D);
ndim           = length(sizeD);

dims = size(D);
if length(dims)==N-1
    dims(end+1) = 1;
end
Rank=size(core);
% [core,Qcore,Pcore] = addRank(core,Qcore,Pcore,Rank);

lambda         = Par.lambda    ;
mu             = Par.mu        ;
rho            = Par.rho       ;
% maxIter        = maxIter-20   ;

maxIter=max(maxIter-10,15);
% disp(maxIter)
%% initialization about M
M         = cell(ndim, 1);
Lam       = cell(ndim, 1);
tempX     = cell(ndim, 1);
sValue    = cell(ndim, 1);
Mnorm     = zeros(ndim, 1);
Msum      = zeros(sizeD);
for i = 1:ndim
    M{i}      = D;
    Lam{i}    = zeros(sizeD);
    tempX{i}  = Unfold(D, sizeD, i);
    sValue{i} = svd(tempX{i}, 'econ');
    Mnorm(i)  = min(sum(sValue{i}>5),Rank(i));
    Msum      = Msum + Fold(M{i},sizeD,i);
end
if length(sizeD)==4
    alpha_c     = circshift(Mnorm, [1,1]).*circshift(Mnorm, [2,2]).*circshift(Mnorm, [3,3]); %computing the  weigth
else
    alpha_c     = circshift(Mnorm, [1,1]).*circshift(Mnorm, [2,2]); %computing the  weigth
end

mu        = mu*(1+1/lambda);
beta      = beta*(1+1/lambda);


%% initialization about other parameters
LamSum    = zeros(sizeD);
temp_n    = zeros(1,ndim);
% Y         = zeros(sizeD);
Y         = D;

%% 
Ps_update=Ps;


%% main loop
for i = 1:floor(maxIter/2.5)
    unfoTemp    = (beta*D+mu*Msum-LamSum)/(beta+ndim*mu);
    %% updating As-N

%     AsT = cellfun(@(x) inv(x.'*x), As(Nlist~=N), 'un', 0);
%     newMultX = tensor(ttm(tensor(unfoTemp), As, -N,'t'));
%     newMultX = tensor(ttm(newMultX, AsT, -N, 't'));
%     newMultX = double(tenmat(newMultX, N));
%     alpha = newMultX*pinv(double(tenmat(core, N)), 1e-8);
    

    newMultX = tensor(ttm(tensor(unfoTemp), As, -N,'t'));
    newMultX = double(tenmat(newMultX, N));
    alpha = newMultX*pinv(double(tenmat(core, N)), 1e-8);
    
    
    
    %% updating As-1:n
    
    for n=1:N-1   
         Xn = reshape(permute(unfoTemp, [n, 1:n-1, n+1:N]), dims(n), []);
         Ps_update{n} = Ps{n}+Xn*kron(alpha, kronC(As(Nlist~=n)))*double(tenmat(core,n)).';
         [V1,~,V2]   = svd(Ps_update{n},'econ');
         As{n}      = V1*V2';
    
    end

    %% updating core

    Qcore_update = Qcore + alpha.'*alpha;
    Pcore_tmp = tensor(ttm(tensor(unfoTemp), [As, {alpha}],'t'));
    Pcore_update = Pcore + Pcore_tmp;
    core = ttm(Pcore_update, inv(Qcore_update), N, 't');
  
    [core,zeroIdxCell]      = ClosedWL1(double(core),1/(beta+ndim*mu),eps);
%     core = tensor(ClosedWL1(double(core),1/(beta+ndim*mu),eps));
    
    %% calculating Y
    Y=double(tensor(ttm(tensor(core), [As, {alpha}])));
    
    %% updating M
    Msum   = 0*Msum;
    LamSum = 0*LamSum;
    for k = 1:ndim
        [tempX{k}, temp_n(k), ~, Mnorm(k)] = Pro2WNNMlogSum(Unfold(Y + Lam{k}/mu, sizeD, k), lambda*alpha_c(k)/mu);%,Rank(k));
        M{k}      = Fold(tempX{k}, sizeD, k);
        Msum      = Msum + M{k};
        if length(sizeD)==4
            alpha_c     = circshift(Mnorm, [1,1]).*circshift(Mnorm, [2,2]).*circshift(Mnorm, [3,3]); %computing the  weigth
        else
            alpha_c     = circshift(Mnorm, [1,1]).*circshift(Mnorm, [2,2]); %computing the  weigth
        end
%         alpha_c     = circshift(Mnorm, [1,1]).*circshift(Mnorm, [2,2]).*circshift(Mnorm, [3,3]); %computing the  weigth
        Lam{k}    = Lam{k}+mu*(Y-M{k});
        LamSum    = LamSum + Lam{k}; % updating the multipliers
    end
%     [core,U,Rank] = ChangerRank(core,U,Rank,zeroIdxCell);
    [core,As,alpha,Pcore,Pcore_update,Qcore,Qcore_update,Ps,Ps_update]=ReduceRk(core,As,alpha,Pcore,Pcore_update,Qcore,Qcore_update,Ps,Ps_update,zeroIdxCell);
    core= tensor(core);
    %% updating mu
    mu         = mu*rho;
    Y_old=Y;
end
end


function [X, n, SigmaNew,wNorm] = Pro2WNNMlogSum(Z, tau)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% min: 1/2*||Z-X||^2 + ||X||_tr
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% [S, V, D, Sigma2] = MySVDtau(Z, tau);
% V = max(diag(V) - tau, 0);
% n = sum(V > 0);
% X = S(:, 1:n) * diag(V(1:n)) * D(:, 1:n)';

%% new
% [m, n] = size(Z);
% if m < n
    AAT = Z*Z';
    [S, Sigma, ~] = svd(AAT);
    Sigma         = sqrt(diag(Sigma));    
    tol           = eps;%max(size(Z)) * eps(max(Sigma));
%     [SigmaNew,n]  = ClosedWNNM(Sigma,tau,tol);
    temp      = (Sigma-tol).^2-4*(tau-tol*Sigma); % log(eps) = 36.0437
    ind    = find (temp>0);
    n         = length(ind);
    % absY      = absY.*ind;
    SigmaNew  = (Sigma(1:n)-tol+sqrt(temp(1:n)))/2;
%     wNorm         = sum(SigmaNew./(Sigma(1:n)+tol));
    wNorm         = sum((log(SigmaNew+eps)+36.0437)/(36.0437+10));
    SigmaNew      = SigmaNew ./ Sigma(1:n);
    X = S(:, 1:n) * diag(SigmaNew) * S(:, 1:n)' * Z;
%     return;
% else
%     [X, n, SigmaNew, wNorm] = Pro2WNNM(Z', tau);
%     X = X';
%     return;
end

function [C,U,UN,Pcore,Pcore_update,Qcore,Qcore_update,Ps,Ps_update] = ReduceRk(C,U,UN,Pcore,Pcore_update,Qcore,Qcore_update,Ps,Ps_update,tempR)

N=numel(U)+1;
nonEmptyIdx = find(~cellfun(@isempty, tempR));
Pcore=double(Pcore);
Pcore_update=double(Pcore_update);
% 遍历非空单元格内容
for i = nonEmptyIdx
    zeroIdx = tempR{i};
%     Rank(i)=Rank(i)-length(zeroIdx);
    if (size(Pcore_update,i)-length(zeroIdx))<2
        NN=size(Pcore_update,i)-2;
        zeroIdx = zeroIdx(1:NN);
    end
    if i<N
        U{i}(:,zeroIdx)=[];
        Ps{i}(:,zeroIdx)=[];
        Ps_update{i}(:,zeroIdx)=[];

    end
    
    switch i
        case 1
            C(zeroIdx, :, :, :) = []; 
            Pcore(zeroIdx, :, :, :) = []; 
            Pcore_update(zeroIdx, :, :, :) = []; 
        case 2
            C(:, zeroIdx, :, :) = []; 
            Pcore(:, zeroIdx, :, :) = []; 
            Pcore_update(:, zeroIdx, :, :) = []; 
        case 3
            if length(size(C))==4
                C(:, :, zeroIdx, :) = []; 
                Pcore(:, :, zeroIdx, :) = []; 
                Pcore_update(:, :, zeroIdx, :) = []; 
            else
                C(:, :,  zeroIdx) = []; 
                Pcore(:, :, zeroIdx) = []; 
                Pcore_update(:, :,  zeroIdx) = []; 
                UN(:, zeroIdx)= []; 
                Qcore(:, zeroIdx)= []; 
                Qcore(zeroIdx,:)= []; 
                Qcore_update(:, zeroIdx)= []; 
                Qcore_update(zeroIdx,:)= []; 
            end
        case 4
            C(:, :, :, zeroIdx) = []; 
            Pcore(:, :, :, zeroIdx) = []; 
            Pcore_update(:, :, :, zeroIdx) = []; 
            UN(:, zeroIdx)= []; 
            Qcore(:, zeroIdx)= []; 
            Qcore(zeroIdx,:)= []; 
            Qcore_update(:, zeroIdx)= []; 
            Qcore_update(zeroIdx,:)= []; 
    end
    
end
Pcore=tensor(Pcore);
end


function [C,Qcore,Pcore] = addRank(C,Qcore,Pcore,newR)
Pcore=double(Pcore);

Rank=size(C);
ind       = true(Rank);
ind_core       = true(Rank(end));

temp      = zeros(newR);
temp(ind) = C;
C         = temp;

temp_Qcore      = zeros(newR(end),newR(end));
temp_Qcore(ind_core) = Qcore;
Qcore=temp_Qcore;


temp_Pcore      = zeros(newR);
temp_Pcore(ind) = Pcore;
Pcore=temp_Pcore;
Pcore=tensor(Pcore);
end
