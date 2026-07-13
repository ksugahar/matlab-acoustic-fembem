function check_projective_geometry_gate(seed)
r=acoustic_fembem.projective_geometry_gate(seed); disp(jsonencode(struct("tool","acoustic_fembem_projective_geometry_gate","ok",r.ok,"result",r)));
if ~r.ok, error("acoustic_fembem:ProjectiveGeometryGateFailed","Projective geometry gate failed."); end
end
