function [mse] = MSE_kron(U,select)
    mse = 1;
    for r=1:length(select)
        T = U{r}(select{r},:)'*U{r}(select{r},:);
        mse = mse * sum(1./eig(T));
    end
end

