function check_text_retrieval_gate(seed)
r=acoustic_fembem.text_retrieval_gate(seed); disp(jsonencode(struct("tool","acoustic_fembem_text_retrieval_gate","ok",r.ok,"result",r)));
if ~r.ok, error("acoustic_fembem:TextRetrievalGateFailed","Text retrieval gate failed."); end
end
