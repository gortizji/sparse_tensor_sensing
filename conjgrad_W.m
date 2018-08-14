function [W] = conjgrad_W(Lw, H, Y, Omega, W)
    W = conjgrad_H(Lw, H, Y', Omega', W);
end