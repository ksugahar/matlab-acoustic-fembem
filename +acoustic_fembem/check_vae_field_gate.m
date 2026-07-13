function check_vae_field_gate(seed)
r=acoustic_fembem.vae_field_gate(seed);
disp(jsonencode(struct("tool","acoustic_fembem_vae_field_gate","ok",r.ok,"result",r)));
if ~r.ok, error("acoustic_fembem:VaeFieldGateFailed","VAE field gate failed."); end
end
