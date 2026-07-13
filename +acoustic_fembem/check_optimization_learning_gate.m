function check_optimization_learning_gate(seed)
r=acoustic_fembem.optimization_learning_gate(seed);
disp(jsonencode(struct("tool","acoustic_fembem_optimization_learning_gate","ok",r.ok,"result",r)));
if ~r.ok, error("acoustic_fembem:OptimizationLearningGateFailed","Optimization learning gate failed."); end
end
