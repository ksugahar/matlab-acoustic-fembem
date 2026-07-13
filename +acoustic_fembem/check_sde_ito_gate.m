function check_sde_ito_gate(seed)
r=acoustic_fembem.sde_ito_gate(seed); disp(jsonencode(struct("tool","acoustic_fembem_sde_ito_gate","ok",r.ok,"result",r)));
if ~r.ok, error("acoustic_fembem:SdeItoGateFailed","SDE Ito gate failed."); end
end
