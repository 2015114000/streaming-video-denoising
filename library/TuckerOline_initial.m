function [ Ps, core, Pcore, Qcore ] = DTuckerO_initial( initX, As, core, R )

% if As is not given, calculate the CP decomposition of the initial data
if ~exist('As') && ~exist('core')
    [cores, factors] = DTucker(tensor(initX), R, 1e-6, 50, 1);
    As = factors;
    core = cores;
    % absorb lambda into the last dimension
end

dims = size(initX);
N = length(dims);

Nlist = 1:N;
for n=1:N-1
%     size(initX)
    Xn = reshape(permute(initX, [n, 1:n-1, n+1:N]), dims(n), []);
%     size(Xn)
    Gn = double(tenmat(core, n));
    Ps{n} = Xn*kronC(As(Nlist~=n))*Gn.';
    AsT = cellfun(@(x,y) x.'*y, As, As, 'un', 0);

end

%%% version 1

AsT = cellfun(@(x) inv(x.'*x), As(Nlist~=N), 'un', 0);
Pcore = tensor(ttm(tensor(initX), As, 't'));
Pcore = tensor(ttm(Pcore, AsT, -N, 't'));
Qcore = As{N}.'*As{N};

end
