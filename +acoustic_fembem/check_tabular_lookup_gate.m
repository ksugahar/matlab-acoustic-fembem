function check_tabular_lookup_gate(seed)
r=acoustic_fembem.tabular_lookup_gate(seed); disp(jsonencode(struct("tool","acoustic_fembem_tabular_lookup_gate","ok",r.ok,"result",r)));
if ~r.ok, error("acoustic_fembem:TabularLookupGateFailed","Tabular lookup gate failed."); end
end
