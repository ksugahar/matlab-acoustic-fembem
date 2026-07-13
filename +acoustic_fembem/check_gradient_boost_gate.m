function check_gradient_boost_gate(seed)
r=acoustic_fembem.gradient_boost_gate(seed); disp(jsonencode(struct("tool","acoustic_fembem_gradient_boost_gate","ok",r.ok,"result",r)));
if ~r.ok, error("acoustic_fembem:GradientBoostGateFailed","Gradient boost gate failed."); end
end
