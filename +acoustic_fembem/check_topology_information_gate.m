function check_topology_information_gate(seed)
r=acoustic_fembem.topology_information_gate(seed); disp(jsonencode(struct("tool","acoustic_fembem_topology_information_gate","ok",r.ok,"result",r)));
if ~r.ok, error("acoustic_fembem:TopologyInformationGateFailed","Topology information gate failed."); end
end
