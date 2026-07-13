function check_policy_trust_gate(seed)
r=acoustic_fembem.policy_trust_gate(seed); disp(jsonencode(struct("tool","acoustic_fembem_policy_trust_gate","ok",r.ok,"result",r)));
if ~r.ok, error("acoustic_fembem:PolicyTrustGateFailed","Policy trust gate failed."); end
end
