function [X,zeroIdxCell]=ClosedWL1(Y,C,oureps)
% solving the following problem
%         sum(w*|Y_i|)+1/2*||Y-X||_F^2
% where w_i =C/(sigmaX_i+oureps),oureps is a small constant
N = length(size(Y));
Nlist = 1:N;


absY      = abs(Y);
signY     = sign(Y);
temp      = (absY-oureps).^2-4*(C-oureps*absY);
ind       = temp>0;
% svp       = sum(ind(:));
% absY      = absY.*ind;

zeroIdxCell=cell(1,N);
for i=1:N
    zeroIdx= find(reshape(all(ind  == 0, Nlist(Nlist~=i)), size(ind,i), 1));
    zeroIdxCell{i}=zeroIdx;
end

absY      = (absY-oureps+sqrt(temp))/2.*ind;
X         = absY.*signY;
end