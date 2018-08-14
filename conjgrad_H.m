function [H] = conjgrad_H(Lh, W, Y, Omega, H)
    b = vec(W'*(Y.*Omega));
    x = vec(H');

    r = b - hessian_mult(Lh,W,x,Omega);
    p = r;
    rsold = r' * r;

    for i = 1:length(b)
        Ap = hessian_mult(Lh,W,p,Omega);
        alpha = rsold / (p' * Ap);
        x = x + alpha * p;
        r = r - alpha * Ap;
        rsnew = r' * r;
        if sqrt(rsnew) < 1e-10
              break;
        end
        p = r + (rsnew / rsold) * p;
        rsold = rsnew;
    end
    H = reshape(x, size(H'))';
end