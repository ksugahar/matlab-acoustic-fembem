function check_graphical_model_gate(seed)
r=acoustic_fembem.graphical_model_gate(seed); disp(jsonencode(struct("tool","acoustic_fembem_graphical_model_gate","ok",r.ok,"result",r)));
if ~r.ok, error("acoustic_fembem:GraphicalModelGateFailed","Graphical model gate failed."); end
end
