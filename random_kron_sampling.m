function [select] = random_kron_sampling(N,b, L)
    R = length(N);
    select = cell(length(N),1);
    DeltaL = L - sum(b);
    DeltaN = N - b;
    
    DeltaLi = randfixedsumint(1,R,DeltaL);
    while any(DeltaLi>DeltaN)
       DeltaLi = randfixedsumint(1,R,DeltaL);
    end
    
    Li = DeltaLi + b;
    for r=1:R
        select{r}=randperm(N(r),Li(r));
    end

end