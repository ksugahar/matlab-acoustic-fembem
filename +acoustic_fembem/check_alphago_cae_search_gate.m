function check_alphago_cae_search_gate(seed)
r=acoustic_fembem.alphago_cae_search_gate(seed);
disp(jsonencode(struct("tool","acoustic_fembem_alphago_cae_search_gate","ok",r.ok,"result",r)));
if ~r.ok, error("acoustic_fembem:AlphaGoCaeSearchGateFailed","AlphaGo CAE search gate failed."); end
end
