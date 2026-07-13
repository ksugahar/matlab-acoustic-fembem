function check_rl_stability_gate(seed)
r=acoustic_fembem.rl_stability_gate(seed);
disp(jsonencode(struct("tool","acoustic_fembem_rl_stability_gate","ok",r.ok,"result",r)));
if ~r.ok, error("acoustic_fembem:RlStabilityGateFailed","RL stability gate failed."); end
end
