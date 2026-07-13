function check_random_matrix_stability_gate(seed)
r=acoustic_fembem.random_matrix_stability_gate(seed); disp(jsonencode(struct("tool","acoustic_fembem_random_matrix_stability_gate","ok",r.ok,"result",r)));
if ~r.ok, error("acoustic_fembem:RandomMatrixStabilityGateFailed","Random matrix stability gate failed."); end
end
