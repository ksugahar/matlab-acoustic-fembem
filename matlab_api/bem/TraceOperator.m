classdef TraceOperator
%TRACEOPERATOR One-hot H1 P1 -> boundary P1 trace injection for a VolMesh.
%
%   T = TraceOperator(mesh);
%   g = T.apply(u);        % or T * u: boundary values of the volume field
%   T.matrix               % sparse (nBem x nFem), exactly one 1 per row
%   T.rowIdentity          % which FEM node and BEM node each row represents
%
% Row k selects volume node femNodeIds(k), so the operator is a pure
% injection: no interpolation, no quadrature. The artifact/schema identity
% strings and the SHA-256 digest of the sparse triplets travel with the
% operator so coupling manifests can bind results to this exact matrix.

properties (Constant)
    operatorPolicy = "one_hot_boundary_node_injection_from_vol_node_ids"
    basisSchemaId = "matlab_h1_p1_to_scalar_bem_p1_trace_basis_v1"
    assemblyRuleId = "first_order_tet_h1_trace_tri_p1_bem_teaching_v1"
    quadratureRuleId = "tri_p1_exact_mass_regular_kernel_teaching_v1"
    observableFamily = "fem_bem_boundary_trace"
end

properties
    matrix        % sparse one-hot trace (nBem x nFem)
    femNodeIds    % volume node id selected by each trace row (nBem x 1)
    bemNodeIds    % compact boundary node ids (1:nBem).'
    rowIdentity   % trace_row_index / fem_node_id / bem_node_id / surface_node_index
    outputDigest  % "sha256:..." of the sparse triplets
    meshId              % identity of the volume mesh
    surfaceMeshId       % identity of the boundary trace mesh
    sourceFileId        % identity of the .vol bytes
    artifactId          % meshId + ":h1_to_scalar_bem_trace_p1"
    operatorArtifactId  % meshId + ":h1_to_scalar_bem_trace_operator_p1"
    outputArtifactId    % meshId + ":h1_to_scalar_bem_trace_output_p1"
    outputPath          % "memory://" + outputArtifactId
    observableId        % meshId + ":h1_to_scalar_bem_boundary_trace_observable_p1"
end

methods
    function T = TraceOperator(mesh)
        arguments
            mesh (1,1) VolMesh
        end
        nFem = size(mesh.vtx, 1);
        nodeIds = mesh.traceNodeIds(:);
        nBem = numel(nodeIds);

        T.matrix = sparse(1:nBem, nodeIds, 1, nBem, nFem);
        T.femNodeIds = nodeIds;
        T.bemNodeIds = (1:nBem).';
        T.rowIdentity = traceRowIdentity(nodeIds);
        T.outputDigest = traceMatrixSha256Id(T.matrix, nBem, nFem);

        T.meshId = mesh.meshId;
        T.surfaceMeshId = mesh.meshId + ":boundary_tri_p1";
        T.sourceFileId = mesh.sourceFileId;
        T.artifactId = mesh.meshId + ":h1_to_scalar_bem_trace_p1";
        T.operatorArtifactId = mesh.meshId + ":h1_to_scalar_bem_trace_operator_p1";
        T.outputArtifactId = mesh.meshId + ":h1_to_scalar_bem_trace_output_p1";
        T.outputPath = "memory://" + T.outputArtifactId;
        T.observableId = mesh.meshId + ":h1_to_scalar_bem_boundary_trace_observable_p1";
    end

    function g = apply(T, u)
        %APPLY Boundary trace of a volume nodal field.
        g = T.matrix * u;
    end

    function g = mtimes(T, u)
        %MTIMES T * u is the same boundary trace as apply(T, u).
        g = T.matrix * u;
    end

    function s = shape(T, dim)
        %SHAPE [nBem, nFem] of the trace matrix (not the object array size).
        s = size(T.matrix);
        if nargin == 2
            s = s(dim);
        end
    end
end
end


function identity = traceRowIdentity(traceNodeIds)
%TRACEROWIDENTITY Bind each trace row to its FEM/BEM node ids.

nRows = numel(traceNodeIds);
identity = struct( ...
    "trace_row_index", num2cell((1:nRows).'), ...
    "fem_node_id", num2cell(traceNodeIds(:)), ...
    "bem_node_id", num2cell(traceNodeIds(:)), ...
    "surface_node_index", num2cell((1:nRows).'));
end


function digest = traceMatrixSha256Id(traceMatrix, nBem, nFem)
%TRACEMATRIXSHA256ID Stable identity for the assembled trace output.

[row, col, val] = find(traceMatrix);
lines = [
    "trace_matrix_p1"
    "n_bem=" + string(nBem)
    "n_fem=" + string(nFem)
    ];
for k = 1:numel(val)
    lines(end + 1, 1) = string(row(k)) + "," + string(col(k)) + "," + sprintf("%.17g", val(k)); %#ok<AGROW>
end
payload = join(lines, newline);
md = javaMethod("getInstance", "java.security.MessageDigest", "SHA-256");
md.update(uint8(char(payload)));
hash = typecast(md.digest(), "uint8");
hex = lower(reshape(dec2hex(hash, 2).', 1, []));
digest = "sha256:" + string(hex);
end
