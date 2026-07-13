function check_pca_covariance_gate(seed)
r=acoustic_fembem.pca_covariance_gate(seed); disp(jsonencode(struct("tool","acoustic_fembem_pca_covariance_gate","ok",r.ok,"result",r)));
if ~r.ok, error("acoustic_fembem:PcaCovarianceGateFailed","PCA covariance gate failed."); end
end
