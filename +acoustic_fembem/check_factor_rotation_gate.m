function check_factor_rotation_gate(seed)
r=acoustic_fembem.factor_rotation_gate(seed); disp(jsonencode(struct("tool","acoustic_fembem_factor_rotation_gate","ok",r.ok,"result",r)));
if ~r.ok, error("acoustic_fembem:FactorRotationGateFailed","Factor rotation gate failed."); end
end
