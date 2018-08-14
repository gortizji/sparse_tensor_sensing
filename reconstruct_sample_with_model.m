function X_est = reconstruct_sample_with_model(X, U1_tilde, U2_tilde, L1, L2)
    Psi1 = U1_tilde(L1, :);
    Psi2 = U2_tilde(L2, :);
    phi_x = X(L2, L1);
    Xf_tilde = pinv(Psi2)*X(L2, L1)*pinv(Psi1).';
    X_est = U2_tilde*Xf_tilde*U1_tilde.';
end