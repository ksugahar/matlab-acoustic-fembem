function check_nice_flow_gate(seed)
r=acoustic_fembem.nice_flow_gate(seed);
disp(jsonencode(struct("tool","acoustic_fembem_nice_flow_gate","ok",r.ok,"result",r)));
if ~r.ok, error("acoustic_fembem:NiceFlowGateFailed","NICE flow gate failed."); end
end
