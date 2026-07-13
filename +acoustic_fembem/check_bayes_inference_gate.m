function check_bayes_inference_gate(seed)
r=acoustic_fembem.bayes_inference_gate(seed); disp(jsonencode(struct("tool","acoustic_fembem_bayes_inference_gate","ok",r.ok,"result",r)));
if ~r.ok, error("acoustic_fembem:BayesInferenceGateFailed","Bayes inference gate failed."); end
end
