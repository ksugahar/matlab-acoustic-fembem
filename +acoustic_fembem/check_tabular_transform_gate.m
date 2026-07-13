function check_tabular_transform_gate(seed)
r=acoustic_fembem.tabular_transform_gate(seed); disp(jsonencode(struct("tool","acoustic_fembem_tabular_transform_gate","ok",r.ok,"result",r)));
if ~r.ok, error("acoustic_fembem:TabularTransformGateFailed","Tabular transform gate failed."); end
end
