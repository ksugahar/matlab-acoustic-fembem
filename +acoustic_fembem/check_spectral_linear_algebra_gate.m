function check_spectral_linear_algebra_gate(seed)
r=acoustic_fembem.spectral_linear_algebra_gate(seed); disp(jsonencode(struct("tool","acoustic_fembem_spectral_linear_algebra_gate","ok",r.ok,"result",r)));
if ~r.ok, error("acoustic_fembem:SpectralLinearAlgebraGateFailed","Spectral linear algebra gate failed."); end
end
