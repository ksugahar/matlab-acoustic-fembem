function fembem_knowledge_tool(topic)
%FEMBEM_KNOWLEDGE_TOOL Print integrated FEM/BEM knowledge for the MCP.
%
% MCP-facing entry point for the MathWorks MATLAB MCP Server custom-tool extension:
% returns one topic of the lab's FEM/BEM knowledge (the integrated Gypsilab
% lane) as plain text so the lab MATLAB MCP knows FEM/BEM the way it knows
% optimization. See acoustic_fembem.fembem_knowledge for the topic list.

arguments
    topic (1,1) string = "overview"
end

disp(acoustic_fembem.fembem_knowledge(topic));
end
