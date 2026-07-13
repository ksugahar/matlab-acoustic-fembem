function check_actor_critic_gate(seed)
r=acoustic_fembem.actor_critic_gate(seed); disp(jsonencode(struct("tool","acoustic_fembem_actor_critic_gate","ok",r.ok,"result",r)));
if ~r.ok, error("acoustic_fembem:ActorCriticGateFailed","Actor critic gate failed."); end
end
