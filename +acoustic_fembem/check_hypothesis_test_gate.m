function check_hypothesis_test_gate(seed)
r=acoustic_fembem.hypothesis_test_gate(seed); disp(jsonencode(struct("tool","acoustic_fembem_hypothesis_test_gate","ok",r.ok,"result",r)));
if ~r.ok, error("acoustic_fembem:HypothesisTestGateFailed","Hypothesis test gate failed."); end
end
