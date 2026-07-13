function check_latent_structural_gate(seed)
r=acoustic_fembem.latent_structural_gate(seed); disp(jsonencode(struct("tool","acoustic_fembem_latent_structural_gate","ok",r.ok,"result",r)));
if ~r.ok, error("acoustic_fembem:LatentStructuralGateFailed","Latent structural gate failed."); end
end
