function check_transformer_scaling_gate(seed)
r=acoustic_fembem.transformer_scaling_gate(seed); disp(jsonencode(struct("tool","acoustic_fembem_transformer_scaling_gate","ok",r.ok,"result",r)));
if ~r.ok, error("acoustic_fembem:TransformerScalingGateFailed","Transformer scaling gate failed."); end
end
