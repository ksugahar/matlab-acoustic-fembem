function check_nlp_transformer_gate(seed)
r=acoustic_fembem.nlp_transformer_gate(seed); disp(jsonencode(struct("tool","acoustic_fembem_nlp_transformer_gate","ok",r.ok,"result",r)));
if ~r.ok, error("acoustic_fembem:NlpTransformerGateFailed","NLP transformer gate failed."); end
end
