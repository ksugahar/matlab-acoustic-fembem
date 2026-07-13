function check_dueling_cae_action_gate(seed)
r=acoustic_fembem.dueling_cae_action_gate(seed);
disp(jsonencode(struct("tool","acoustic_fembem_dueling_cae_action_gate","ok",r.ok,"result",r)));
if ~r.ok, error("acoustic_fembem:DuelingCaeActionGateFailed","Dueling CAE action gate failed."); end
end
