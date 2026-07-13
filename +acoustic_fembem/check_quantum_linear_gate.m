function check_quantum_linear_gate(seed)
r=acoustic_fembem.quantum_linear_gate(seed); disp(jsonencode(struct("tool","acoustic_fembem_quantum_linear_gate","ok",r.ok,"result",r)));
if ~r.ok, error("acoustic_fembem:QuantumLinearGateFailed","Quantum linear gate failed."); end
end
