function X=my_normalized(X)
X=double(X);
Nway = size(X);
ndims=length(Nway);
if ndims==3
    for i=1:Nway(3)
        X(:,:,i) = X(:,:,i)/max(max(X(:,:,i)));
    end
else
    for i=1:Nway(4)
        X(:,:,:,i) = X(:,:,:,i)/max(max(max(X(:,:,:,i))));

    end
end


% Nway = size(X);
% for i=1:Nway(3)
%     X(:,:,i) = (X(:,:,i)-min(min(X(:,:,i))))/(max(max(X(:,:,i)))-min(min(X(:,:,i))));
% end
% end