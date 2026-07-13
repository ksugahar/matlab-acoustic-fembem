function check_latent_distribution_gate(seed)
r=acoustic_fembem.latent_distribution_gate(seed); disp(jsonencode(struct("tool","acoustic_fembem_latent_distribution_gate","ok",r.ok,"result",r)));
if ~r.ok, error("acoustic_fembem:LatentDistributionGateFailed","Latent distribution gate failed."); end
end
