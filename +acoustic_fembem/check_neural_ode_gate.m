function check_neural_ode_gate(seed)
r=acoustic_fembem.neural_ode_gate(seed); disp(jsonencode(struct("tool","acoustic_fembem_neural_ode_gate","ok",r.ok,"result",r)));
if ~r.ok, error("acoustic_fembem:NeuralOdeGateFailed","Neural ODE gate failed."); end
end
