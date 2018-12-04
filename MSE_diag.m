function [mse] = MSE_diag(U,select)
    T = eye(size(U,2), size(U,2));
    for r=1:length(select)
        T = T .* (U{r}(select{r},:)'*U{r}(select{r},:));
    end
    mse = sum(1./eig(T));
end

