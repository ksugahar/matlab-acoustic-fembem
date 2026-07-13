function check_activation_probability_gate(seed)
r=acoustic_fembem.activation_probability_gate(seed); disp(jsonencode(struct("tool","acoustic_fembem_activation_probability_gate","ok",r.ok,"result",r)));
if ~r.ok, error("acoustic_fembem:ActivationProbabilityGateFailed","Activation probability gate failed."); end
end
