function check_diffusion_inverse_gate(seed)
r=acoustic_fembem.diffusion_inverse_gate(seed);
disp(jsonencode(struct("tool","acoustic_fembem_diffusion_inverse_gate","ok",r.ok,"result",r)));
if ~r.ok, error("acoustic_fembem:DiffusionInverseGateFailed","Diffusion inverse gate failed."); end
end
