function check_time_series_gate(seed)
r=acoustic_fembem.time_series_gate(seed); disp(jsonencode(struct("tool","acoustic_fembem_time_series_gate","ok",r.ok,"result",r)));
if ~r.ok, error("acoustic_fembem:TimeSeriesGateFailed","Time series gate failed."); end
end
