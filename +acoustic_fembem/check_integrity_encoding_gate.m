function check_integrity_encoding_gate(seed)
r=acoustic_fembem.integrity_encoding_gate(seed); disp(jsonencode(struct("tool","acoustic_fembem_integrity_encoding_gate","ok",r.ok,"result",r)));
if ~r.ok, error("acoustic_fembem:IntegrityEncodingGateFailed","Integrity encoding gate failed."); end
end
