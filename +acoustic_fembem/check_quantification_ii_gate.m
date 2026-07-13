function check_quantification_ii_gate(seed)
r=acoustic_fembem.quantification_ii_gate(seed); disp(jsonencode(struct("tool","acoustic_fembem_quantification_ii_gate","ok",r.ok,"result",r)));
if ~r.ok, error("acoustic_fembem:QuantificationIiGateFailed","Quantification II gate failed."); end
end
