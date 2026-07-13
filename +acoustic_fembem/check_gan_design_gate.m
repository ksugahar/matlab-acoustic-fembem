function check_gan_design_gate(seed)
r=acoustic_fembem.gan_design_gate(seed);
disp(jsonencode(struct("tool","acoustic_fembem_gan_design_gate","ok",r.ok,"result",r)));
if ~r.ok, error("acoustic_fembem:GanDesignGateFailed","GAN design gate failed."); end
end
