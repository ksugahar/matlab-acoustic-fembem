function check_recurrent_sequence_gate(seed)
r=acoustic_fembem.recurrent_sequence_gate(seed); disp(jsonencode(struct("tool","acoustic_fembem_recurrent_sequence_gate","ok",r.ok,"result",r)));
if ~r.ok, error("acoustic_fembem:RecurrentSequenceGateFailed","Recurrent sequence gate failed."); end
end
