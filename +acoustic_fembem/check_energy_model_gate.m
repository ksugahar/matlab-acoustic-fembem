function check_energy_model_gate(seed)
r=acoustic_fembem.energy_model_gate(seed);
disp(jsonencode(struct("tool","acoustic_fembem_energy_model_gate","ok",r.ok,"result",r)));
if ~r.ok, error("acoustic_fembem:EnergyModelGateFailed","Energy model gate failed."); end
end
