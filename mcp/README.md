# Acoustic FEM/BEM MCP Layer

This directory exposes the MATLAB acoustic FEM/BEM teaching code through the
official MathWorks MATLAB MCP Server.

The solver, validation fixtures, and MCP entry points live in the same
repository:

- MATLAB package: `+acoustic_fembem`
- Extension file: `mcp/extensions/acoustic-fembem-tools.json`
- Verifier: `mcp/tools/verify_acoustic_fembem_mcp.ps1`

## Tools

| Tool | MATLAB function | Purpose |
|------|-----------------|---------|
| `acoustic_fembem_activation_probability_gate` | `acoustic_fembem.check_activation_probability_gate` | Sigmoid/softmax gradient and saturation checks |
| `acoustic_fembem_topology_information_gate` | `acoustic_fembem.check_topology_information_gate` | Tropical/TDA/manifold/Fisher metric checks |
| `acoustic_fembem_integrity_encoding_gate` | `acoustic_fembem.check_integrity_encoding_gate` | Chained integrity and UTF-8 binary controls |
| `acoustic_fembem_tabular_transform_gate` | `acoustic_fembem.check_tabular_transform_gate` | String extraction and conditional aggregation checks |
| `acoustic_fembem_bayes_inference_gate` | `acoustic_fembem.check_bayes_inference_gate` | Conjugacy, MLE/MAP, MCMC, and R-squared controls |
| `acoustic_fembem_time_series_gate` | `acoustic_fembem.check_time_series_gate` | VAR, impulse, Granger, and innovation checks |
| `acoustic_fembem_gradient_boost_gate` | `acoustic_fembem.check_gradient_boost_gate` | Residual boosting and overstep negative control |
| `acoustic_fembem_neural_ode_gate` | `acoustic_fembem.check_neural_ode_gate` | ODE convergence and sensitivity checks |
| `acoustic_fembem_tabular_lookup_gate` | `acoustic_fembem.check_tabular_lookup_gate` | Exact/approximate lookup and missing-key controls |
| `acoustic_fembem_graphical_model_gate` | `acoustic_fembem.check_graphical_model_gate` | Conditional independence, precision, and DAG controls |
| `acoustic_fembem_sde_ito_gate` | `acoustic_fembem.check_sde_ito_gate` | Ito Monte Carlo/isometry and Black-Scholes PDE checks |
| `acoustic_fembem_cnn_architecture_gate` | `acoustic_fembem.check_cnn_architecture_gate` | CNN shape, residual, separable, Inception and SENet checks |
| `acoustic_fembem_recurrent_sequence_gate` | `acoustic_fembem.check_recurrent_sequence_gate` | RNN/LSTM/GRU sequence and stability checks |
| `acoustic_fembem_latent_structural_gate` | `acoustic_fembem.check_latent_structural_gate` | Hierarchical Bayes, SEM, and categorical response checks |
| `acoustic_fembem_quantification_ii_gate` | `acoustic_fembem.check_quantification_ii_gate` | Generalized eigenvalue and discrimination-ratio checks |
| `acoustic_fembem_factor_rotation_gate` | `acoustic_fembem.check_factor_rotation_gate` | Orthogonal/oblique factor-rotation checks |
| `acoustic_fembem_text_retrieval_gate` | `acoustic_fembem.check_text_retrieval_gate` | Embedding, tf-idf, and BM25 ranking controls |
| `acoustic_fembem_projective_geometry_gate` | `acoustic_fembem.check_projective_geometry_gate` | Homogeneous and conic-invariance checks |
| `acoustic_fembem_nlp_transformer_gate` | `acoustic_fembem.check_nlp_transformer_gate` | Attention, masks, sharing, BLEU, and text-task schema |
| `acoustic_fembem_quantification_gate` | `acoustic_fembem.check_quantification_gate` | Preference SVD and constrained eigenvalue checks |
| `acoustic_fembem_hypothesis_test_gate` | `acoustic_fembem.check_hypothesis_test_gate` | Null/alternative and p-value negative controls |
| `acoustic_fembem_transformer_scaling_gate` | `acoustic_fembem.check_transformer_scaling_gate` | Scaling law, sparse causality, and extrapolation withholding |
| `acoustic_fembem_pca_covariance_gate` | `acoustic_fembem.check_pca_covariance_gate` | PCA variance/reconstruction and covariance geometry |
| `acoustic_fembem_quantum_linear_gate` | `acoustic_fembem.check_quantum_linear_gate` | Unitarity, Deutsch control, and I/O-cost checks |
| `acoustic_fembem_bellman_control_gate` | `acoustic_fembem.check_bellman_control_gate` | Bellman, TD, SARSA/Q-learning, and GPI checks |
| `acoustic_fembem_spectral_linear_algebra_gate` | `acoustic_fembem.check_spectral_linear_algebra_gate` | SVD, isometry, diagonalization, and covariance checks |
| `acoustic_fembem_random_matrix_stability_gate` | `acoustic_fembem.check_random_matrix_stability_gate` | Deep-Jacobian vanishing/explosion negative controls |
| `acoustic_fembem_latent_distribution_gate` | `acoustic_fembem.check_latent_distribution_gate` | Generation/identification probability and solver checks |
| `acoustic_fembem_linear_inverse_gate` | `acoustic_fembem.check_linear_inverse_gate` | Pseudoinverse, regression, TSVD, and CCA checks |
| `acoustic_fembem_actor_critic_gate` | `acoustic_fembem.check_actor_critic_gate` | REINFORCE, actor-critic, trace, and DQN target checks |
| `acoustic_fembem_a3c_async_gate` | `acoustic_fembem.check_a3c_async_gate` | Async worker version/staleness, on-policy, and solver-promotion checks |
| `acoustic_fembem_policy_trust_gate` | `acoustic_fembem.check_policy_trust_gate` | DPG gradient, TRPO KL limit, and DDPG target controls |
| `acoustic_fembem_variational_geometry_gate` | `acoustic_fembem.check_variational_geometry_gate` | KL/ELBO/free-energy and geometry-preservation checks |
| `acoustic_fembem_rl_stability_gate` | `acoustic_fembem.check_rl_stability_gate` | Prioritized replay correction and Double-DQN bias controls |
| `acoustic_fembem_energy_model_gate` | `acoustic_fembem.check_energy_model_gate` | Exact Boltzmann/RBM gradient and DBN handoff checks |
| `acoustic_fembem_dueling_cae_action_gate` | `acoustic_fembem.check_dueling_cae_action_gate` | Centered Q=V+A, action contract, and solver reward replay |
| `acoustic_fembem_vae_field_gate` | `acoustic_fembem.check_vae_field_gate` | ELBO/KL, reparameterization gradient, OOD, and forward-QoI gate |
| `acoustic_fembem_gan_design_gate` | `acoustic_fembem.check_gan_design_gate` | Mode coverage, collapse rejection, discriminator, and solver promotion |
| `acoustic_fembem_alphago_cae_search_gate` | `acoustic_fembem.check_alphago_cae_search_gate` | Policy/value/tree-search separation with solver verification |
| `acoustic_fembem_nice_flow_gate` | `acoustic_fembem.check_nice_flow_gate` | Reversible coupling, Jacobian, likelihood, OOD, and forward-QoI replay |
| `acoustic_fembem_diffusion_inverse_gate` | `acoustic_fembem.check_diffusion_inverse_gate` | Diffusion noise schedule, least-squares denoising, and forward-residual gate |
| `acoustic_fembem_optimization_learning_gate` | `acoustic_fembem.check_optimization_learning_gate` | Gradient, GP/Bayes, rank selection, PSO/CMA schema, and forward replay gate |
| `acoustic_fembem_ml_cae_pipeline` | `acoustic_fembem.check_ml_cae_pipeline` | Reproducible schema, POD, GP uncertainty, active learning, solver verification, and generative/RL preflight |
| `acoustic_fembem_check_result_manifest_file` | `acoustic_fembem.check_result_manifest_file` | Validate script/MCP-ready result manifests |
| `acoustic_fembem_acoustic_gate` | `acoustic_fembem.check_fembem_acoustic_gate` | Run acoustic FEM/BEM gates against analytic references |
| `acoustic_fembem_crossval_gate` | `acoustic_fembem.check_fembem_crossval_gate` | Run `.vol`-backed cross-validation against radia-ngsolve/NGSolve references |
| `acoustic_fembem_knowledge` | `acoustic_fembem.fembem_knowledge_tool` | Serve compact acoustic FEM/BEM teaching knowledge |
| `acoustic_fembem_vol_mesh_summary` | `acoustic_fembem.check_vol_mesh_summary` | Summarize `.vol` meshes and viewer guidance |
| `acoustic_fembem_hmatrix_scaling` | `acoustic_fembem.check_hmatrix_scaling` | Run readable ACA+ scaling and gate dense-reference error, rank, and storage growth |
| `acoustic_fembem_repository_health` | `acoustic_fembem.check_repository_health` | Pre-push health check for the integrated repo and MCP extension |
| `acoustic_fembem_rwg_hcurl_trace_artifact_gate` | `acoustic_fembem_check_rwg_hcurl_trace_artifact_gate` | Gate GYP-081..090 RWG/HCurl trace, de Rham, and independent reference evidence |

The acoustic gate includes adjoint optimization cases:

- `focus_adjoint`: gradient of focused pressure intensity vs finite difference.
- `radiation_force`: acoustic radiation-force postprocess from the BEM field.
- `thrust_adjoint`: elastic-bead force quadratic form and Wirtinger gradient.

The knowledge tool includes `matlab_execution_policy`, `vol_visualization`, and
`pde_vol_bridge`.  The policy is normal `.m` functions/scripts plus MCP JSON,
Netgen for native interactive `.vol` inspection, MATLAB `plotVolMesh` for
figure previews, and `acoustic_fembem_vol_mesh_summary` for LLM/headless
preflight.

## Official Server

Install the official MathWorks MATLAB MCP Server separately, then pass:

```powershell
--extension-file=<repo>\mcp\extensions\acoustic-fembem-tools.json
```

For an already shared R2026a session, preflight that exactly one shared
session is visible, then use the official server's capture route:

```powershell
--matlab-session-mode=auto --matlab-display-mode=nodesktop `
  --extension-file=<repo>\mcp\extensions\acoustic-fembem-tools.json
```

The official CLI does not allow `existing` together with `nodesktop`.  The
`auto` route is acceptable here only after the single-session preflight.  The
MCP log must say `Attaching to existing session`, and the MATLAB process set
must be unchanged before and after the call.  This keeps command-window output
on the capture `FEval` path without closing or replacing the shared desktop
session.

The official server supplies the MCP runtime.  This repository supplies the
domain tools and validation gates.

## MATLAB Agentic Toolkit

Use the MathWorks MATLAB Agentic Toolkit as the setup and skills reference
layer.  It can install/update the official MATLAB MCP Server and register
agent skills for MATLAB workflows.  This repository does not replace those
pieces; it adds only the acoustic FEM/BEM/CQ domain extension.

Recommended split:

- official MATLAB MCP Server: MATLAB session/runtime, code execution, tests,
  Code Analyzer checks, and toolbox detection;
- MATLAB Agentic Toolkit: skills and setup guidance, selected narrowly for the
  project;
- `matlab-acoustic-fembem`: `.vol` tri/tet intake, P1 FEM/BEM, Johnson-Nedelec
  coupling, Lubich CQ, Gmsh artifacts, and NGSolve.BEM/radia cross-validation
  gates.

When an external solver bridge is involved, attach to the already shared MATLAB
session through the official server's existing-session workflow. Pure MATLAB
acoustic FEM/BEM checks may run separately when they do not touch that bridge.

## Repository Boundary

The acoustic FEM/BEM MCP layer stays inside this repository.  Do not maintain a
separate lab `matlab-mcp` repository for these tools: the extension file is the
public contract, while the official MathWorks server remains the runtime.
Internal cross-validation provenance, local solver logs, and private MATLAB
automation notes stay outside the public package unless rewritten as scrubbed
domain knowledge.
