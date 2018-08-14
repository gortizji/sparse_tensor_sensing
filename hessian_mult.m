function [v] = hessian_mult(Lh, W, s, Omega)
    S = reshape(s, size(W,2), []);
    K = zeros(size(S));
    
    WS_omega = (W*S).*Omega;
    K = W' * WS_omega;
    V = K + S*Lh;
    v = vec(V);
end