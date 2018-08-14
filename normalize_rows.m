function [U] = normalize_rows(U)

for i=1:size(U,1)
    U(i,:)=U(i,:)/norm(U(i,:),2);
end
end

