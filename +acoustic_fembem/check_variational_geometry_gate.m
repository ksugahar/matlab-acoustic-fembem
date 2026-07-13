function check_variational_geometry_gate(seed)
r=acoustic_fembem.variational_geometry_gate(seed); disp(jsonencode(struct("tool","acoustic_fembem_variational_geometry_gate","ok",r.ok,"result",r)));
if ~r.ok, error("acoustic_fembem:VariationalGeometryGateFailed","Variational geometry gate failed."); end
end
