function resultMat = kronC(factors)

N = length(factors);

if N == 1
    
   resultMat = factors{1};
    
else

    tmp = kron(factors{N}, factors{N-1});
    for n=N-2:-1:1
        tmp = kron(tmp, factors{n});
    end

    resultMat = tmp;

end

