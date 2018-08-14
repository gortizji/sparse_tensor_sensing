function [Y, W, H] = grals(Lh, Lw, Omega, Y, k, Otest)
%GRALS Graph Regularized Alternating Least Squares
% Implementation by Guillermo Ortiz Jimenez (g.ortizjimenez@student.tudelft.nl)
% Algorithm based on paper:
% Rao, N., Yu, H. F., Ravikumar, P. K., & Dhillon, I. S. (2015). 
% "Collaborative filtering with graph information: Consistency 
% and scalable methods". In Advances in neural information processing systems 
% (pp. 2107-2115).
% Inputs:
%   Lh: Columns laplacian
%   Lw: Rows laplacian
%   Omega: Observed training data mask
%   Y: Measured data
%   k: Rank of decomposition
%   Otest: Observed test mask
% Outputs:
%   Y: Estimated data (Y=W*H')
%   W: Row factor matrix
%   H: Column factor matrix

    m = size(Omega,1);
    n = size(Omega,2);
    
    [W, ~, H] = svd(Y);
    W = W(:,1:k);
    H = H(:,1:k);
    Y_ = W*H';
    for n = 1:10
        H = conjgrad_H(Lh, W, Y, Omega, H);
        W = conjgrad_W(Lw, H, Y, Omega, W);
        Y_ = W*H';
        e = sqrt(norm(Otest.*Y - Otest.*Y_,'fro')^2/sum(Otest(:)))
    end
    Y = W*H';
end

