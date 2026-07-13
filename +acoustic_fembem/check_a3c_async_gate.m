function check_a3c_async_gate(seed)
r=acoustic_fembem.a3c_async_gate(seed);
disp(jsonencode(struct("tool","acoustic_fembem_a3c_async_gate","ok",r.ok,"result",r)));
end
