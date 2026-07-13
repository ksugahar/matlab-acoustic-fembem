function check_linear_inverse_gate(seed)
r=acoustic_fembem.linear_inverse_gate(seed); disp(jsonencode(struct("tool","acoustic_fembem_linear_inverse_gate","ok",r.ok,"result",r)));
if ~r.ok, error("acoustic_fembem:LinearInverseGateFailed","Linear inverse gate failed."); end
end
