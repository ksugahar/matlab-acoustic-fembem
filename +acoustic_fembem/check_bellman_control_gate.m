function check_bellman_control_gate(seed)
r=acoustic_fembem.bellman_control_gate(seed); disp(jsonencode(struct("tool","acoustic_fembem_bellman_control_gate","ok",r.ok,"result",r)));
if ~r.ok, error("acoustic_fembem:BellmanControlGateFailed","Bellman control gate failed."); end
end
