function check_cnn_architecture_gate(seed)
r=acoustic_fembem.cnn_architecture_gate(seed); disp(jsonencode(struct("tool","acoustic_fembem_cnn_architecture_gate","ok",r.ok,"result",r)));
if ~r.ok, error("acoustic_fembem:CnnArchitectureGateFailed","CNN architecture gate failed."); end
end
