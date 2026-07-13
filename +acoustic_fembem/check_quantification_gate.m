function check_quantification_gate(seed)
r=acoustic_fembem.quantification_gate(seed); disp(jsonencode(struct("tool","acoustic_fembem_quantification_gate","ok",r.ok,"result",r)));
if ~r.ok, error("acoustic_fembem:QuantificationGateFailed","Quantification gate failed."); end
end
