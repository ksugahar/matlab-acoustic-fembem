function tests = testMcpAcousticFembemTools
%TESTMCPACOUSTICFEMBEMTOOLS MCP-facing acoustic FEM/BEM entry points.

tests = functiontests(localfunctions);
end


function setupOnce(testCase)
repoRoot = fileparts(fileparts(mfilename("fullpath")));
testCase.applyFixture(matlab.unittest.fixtures.PathFixture(repoRoot));
addpath(genpath(fullfile(repoRoot, "matlab_api")));
addpath(fullfile(repoRoot, "examples"));
addpath(fullfile(repoRoot, "validation"));
end


function testRepositoryRoot(testCase)
root = acoustic_fembem.repository_root();
verifyTrue(testCase, isfolder(fullfile(root, "matlab_api")));
verifyTrue(testCase, isfile(fullfile(root, "mcp", "extensions", "acoustic-fembem-tools.json")));
end


function testMcpLayerStaysInsideAcousticFembemRepo(testCase)
root = acoustic_fembem.repository_root();
readme = string(fileread(fullfile(root, "mcp", "README.md")));
requirements = string(fileread(fullfile(root, "mcp", "REQUIREMENTS.md")));
verifySubstring(testCase, readme, "stays inside this repository");
verifySubstring(testCase, readme, "official MathWorks server remains the runtime");
verifySubstring(testCase, readme, "MATLAB Agentic Toolkit");
verifySubstring(testCase, readme, "existing-session workflow");
verifySubstring(testCase, requirements, "extension is intentionally thin");
verifySubstring(testCase, requirements, "process/session management");
verifySubstring(testCase, requirements, "matlab-agentic-toolkit");
verifySubstring(testCase, requirements, "Existing MATLAB Session Policy");
end


function testKnowledgeIncludesCrossvalTopic(testCase)
body = acoustic_fembem.fembem_knowledge("radia_ngsolve_crossval");
verifyGreaterThan(testCase, strlength(body), 300);
verifySubstring(testCase, body, ".vol");
verifySubstring(testCase, body, "radia-ngsolve");
end


function testKnowledgeIncludesNgsolveBem50Topic(testCase)
body = acoustic_fembem.fembem_knowledge("ngsolve_bem_50");
verifyGreaterThan(testCase, strlength(body), 600);
verifySubstring(testCase, body, "50 examples");
verifySubstring(testCase, body, "NGSolve.BEM");
verifySubstring(testCase, body, "Netgen");
verifySubstring(testCase, body, "Do not claim an analytic solution");
end


function testKnowledgeIncludesPdeVolBridgeTopic(testCase)
body = acoustic_fembem.fembem_knowledge("pde_vol_bridge");
verifyGreaterThan(testCase, strlength(body), 300);
verifySubstring(testCase, body, "PDE Toolbox");
verifySubstring(testCase, body, "writePdeMeshVol");
verifySubstring(testCase, body, ".vol");
end


function testKnowledgeIncludesVolVisualizationTopic(testCase)
body = acoustic_fembem.fembem_knowledge("vol_visualization");
verifyGreaterThan(testCase, strlength(body), 300);
verifySubstring(testCase, body, "Netgen");
verifySubstring(testCase, body, "plotVolMesh");
verifySubstring(testCase, body, "acoustic_fembem_vol_mesh_summary");
verifySubstring(testCase, body, "Ng_LoadMesh");
verifySubstring(testCase, body, ".sol files are mesh-free");
end


function testKnowledgeIncludesConvolutionQuadratureTopic(testCase)
body = acoustic_fembem.fembem_knowledge("convolution_quadrature");
verifyGreaterThan(testCase, strlength(body), 400);
verifySubstring(testCase, body, "Lubich CQ");
verifySubstring(testCase, body, "laplaceSingleLayerGalerkin");
verifySubstring(testCase, body, "A-stability");
verifySubstring(testCase, body, "IMAGINARY axis");
verifyEqual(testCase, acoustic_fembem.fembem_knowledge("cq"), body);
end


function testCqTimeGridMcpContract(testCase)
report = cqTimeGridManifest(1e-3, 100, "BDF2", 128, 0);
verifyTrue(testCase, report.ok);
verifyEqual(testCase, report.contourSamples, 128);
verifyEqual(testCase, report.timeEnd, 0.099, "AbsTol", 1e-14);
verifyGreaterThan(testCase, report.minRealLaplaceNode, 0);

bad = cqTimeGridManifest(1e-3, 100, "BDF2", 64, 0.9);
verifyFalse(testCase, bad.ok);
verifyFalse(testCase, bad.checks.contourCoversTimeGrid);

out = evalc("acoustic_fembem.check_cq_time_grid(1e-3, 100, ""BDF2"", 128, 0)");
decoded = jsondecode(out);
verifyTrue(testCase, decoded.ok);
verifyEqual(testCase, string(decoded.tool), "acoustic_fembem_cq_time_grid");
verifyEqual(testCase, string(decoded.result.schema), "matlab-acoustic-fembem.cq-time-grid.v1");
end


function testBalancedLearningProfileMcpContract(testCase)
profile = balancedLearningProfile();
verifyEqual(testCase, profile.policy, "equal_capability_gain_v1");
verifyEqual(testCase, profile.stage_count, 10);
verifyEqual(testCase, numel(unique([profile.stages.capability_id])), 10);
verifyEqual(testCase, sort(string(fieldnames(profile.workflow_roles))), ...
    sort(["detect"; "check"; "run"; "test"]));
verifyEqual(testCase, profile.self_check.status, "ok");

bad = profile;
bad.stages(4).negative_control = "";
rejected = validateBalancedLearningProfile(bad);
verifyEqual(testCase, rejected.status, "needs_attention");
verifyFalse(testCase, rejected.checks.controls_complete);

out = evalc("acoustic_fembem.check_balanced_learning_profile()");
decoded = jsondecode(out);
verifyTrue(testCase, decoded.ok);
verifyEqual(testCase, string(decoded.tool), "acoustic_fembem_balanced_learning_profile");
end


function testMlCaePipelineContract(testCase)
report = acoustic_fembem.ml_cae_pipeline("all", 17);
verifyTrue(testCase, report.ok);
verifyEqual(testCase, report.seed, 17);
verifyFalse(testCase, report.provenance.generated_candidate_is_ground_truth);
verifyLessThan(testCase, report.pod.orthogonality_error, 1e-10);
verifyGreaterThanOrEqual(testCase, report.active_learning.next_x, 0);
verifyLessThanOrEqual(testCase, report.active_learning.next_x, 1);
verifyTrue(testCase, report.forward_verification.passed);
verifyTrue(testCase, report.generative_preflight.requires_forward_verification);
verifyTrue(testCase, report.rl_preflight.requires_forward_verification);

out = evalc("acoustic_fembem.check_ml_cae_pipeline(""forward_verify"", 17)");
decoded = jsondecode(out);
verifyTrue(testCase, decoded.ok);
verifyEqual(testCase, string(decoded.tool), "acoustic_fembem_ml_cae_pipeline");
verifyTrue(testCase, decoded.result.forward_verification.candidate_promoted);
end


function testMlCaePipelineRejectsUnknownStage(testCase)
verifyError(testCase, @() acoustic_fembem.ml_cae_pipeline("unknown", 17), ...
    "acoustic_fembem:UnknownMlCaeStage");
end


function testOptimizationLearningGate(testCase)
r = acoustic_fembem.optimization_learning_gate(23);
verifyTrue(testCase,r.ok);
verifyLessThan(testCase,r.gradient.complex_step_relative_error,1e-12);
verifyTrue(testCase,r.pso.bounds_ok);
verifyTrue(testCase,r.cma_es.bounds_ok);
verifyEqual(testCase,r.pso.history_schema,r.cma_es.history_schema);
verifyTrue(testCase,r.forward_verification.passed);
out=evalc("acoustic_fembem.check_optimization_learning_gate(23)");
j=jsondecode(out); verifyTrue(testCase,j.ok);
verifyEqual(testCase,string(j.tool),"acoustic_fembem_optimization_learning_gate");
end


function testDiffusionInverseGate(testCase)
a=acoustic_fembem.diffusion_inverse_gate(31);
b=acoustic_fembem.diffusion_inverse_gate(31);
verifyTrue(testCase,a.ok);
verifyEqual(testCase,a.rows,b.rows,"AbsTol",0);
verifyFalse(testCase,a.generated_candidate_is_ground_truth);
verifyTrue(testCase,a.promotion_requires_forward_solver);
verifyLessThan(testCase,mean(a.rows(:,4)),mean(a.rows(:,3)));
out=evalc("acoustic_fembem.check_diffusion_inverse_gate(31)");
j=jsondecode(out); verifyTrue(testCase,j.ok);
verifyEqual(testCase,string(j.tool),"acoustic_fembem_diffusion_inverse_gate");
end


function testNiceFlowGate(testCase)
r=acoustic_fembem.nice_flow_gate(37);
verifyTrue(testCase,r.ok);
verifyLessThan(testCase,r.roundtrip_max_abs,1e-12);
verifyEqual(testCase,r.jacobian_determinant,1,"AbsTol",1e-12);
verifyLessThan(testCase,r.mean_nll_after,r.mean_nll_before);
verifyFalse(testCase,r.generated_candidate_is_ground_truth);
out=evalc("acoustic_fembem.check_nice_flow_gate(37)"); j=jsondecode(out);
verifyTrue(testCase,j.ok); verifyEqual(testCase,string(j.tool),"acoustic_fembem_nice_flow_gate");
end


function testAlphaGoCaeSearchGate(testCase)
r=acoustic_fembem.alphago_cae_search_gate(41);
verifyTrue(testCase,r.ok); verifyEqual(testCase,sum(r.visit_counts_root),r.simulations);
verifyEqual(testCase,r.candidate,r.global_design); verifyFalse(testCase,r.candidate_is_ground_truth);
out=evalc("acoustic_fembem.check_alphago_cae_search_gate(41)"); j=jsondecode(out);
verifyTrue(testCase,j.ok); verifyEqual(testCase,string(j.tool),"acoustic_fembem_alphago_cae_search_gate");
end


function testGanDesignGate(testCase)
r=acoustic_fembem.gan_design_gate(43);
verifyTrue(testCase,r.ok); verifyEqual(testCase,r.good_mode_coverage,r.expected_modes);
verifyEqual(testCase,r.collapsed_mode_coverage,1);
verifyGreaterThan(testCase,r.collapsed_discriminator_accuracy,r.good_discriminator_accuracy);
verifyFalse(testCase,r.generated_candidate_is_ground_truth);
out=evalc("acoustic_fembem.check_gan_design_gate(43)"); j=jsondecode(out);
verifyTrue(testCase,j.ok); verifyEqual(testCase,string(j.tool),"acoustic_fembem_gan_design_gate");
end


function testVaeFieldGate(testCase)
r=acoustic_fembem.vae_field_gate(47);
verifyTrue(testCase,r.ok); verifyGreaterThanOrEqual(testCase,r.kl_mean,0);
verifyLessThan(testCase,r.reconstruction_mse,r.mean_baseline_mse);
verifyLessThan(testCase,r.forward_qoi_relative_error,r.baseline_qoi_relative_error);
verifyLessThan(testCase,r.reparameterization_gradient_error,1e-6);
verifyFalse(testCase,r.generated_candidate_is_ground_truth);
out=evalc("acoustic_fembem.check_vae_field_gate(47)"); j=jsondecode(out);
verifyTrue(testCase,j.ok); verifyEqual(testCase,string(j.tool),"acoustic_fembem_vae_field_gate");
end


function testDuelingCaeActionGate(testCase)
r=acoustic_fembem.dueling_cae_action_gate(53);
verifyTrue(testCase,r.ok);
verifyLessThan(testCase,max(abs(mean(r.advantage,2))),1e-12);
verifyTrue(testCase,r.checks.offset_invariant_advantage);
verifyTrue(testCase,r.checks.action_contract_rejects_invalid);
verifyFalse(testCase,r.policy_candidate_is_ground_truth);
out=evalc("acoustic_fembem.check_dueling_cae_action_gate(53)"); j=jsondecode(out);
verifyTrue(testCase,j.ok); verifyEqual(testCase,string(j.tool),"acoustic_fembem_dueling_cae_action_gate");
end


function testRlStabilityGate(testCase)
r=acoustic_fembem.rl_stability_gate(59); verifyTrue(testCase,r.ok);
verifyGreaterThan(testCase,r.dqn_overestimation_bias,abs(r.double_dqn_bias));
verifyEqual(testCase,sum(r.sampling_probability),1,"AbsTol",1e-12);
out=evalc("acoustic_fembem.check_rl_stability_gate(59)"); j=jsondecode(out); verifyTrue(testCase,j.ok);
end


function testEnergyModelGate(testCase)
r=acoustic_fembem.energy_model_gate(61); verifyTrue(testCase,r.ok);
verifyGreaterThan(testCase,r.log_likelihood_after,r.log_likelihood_before);
verifyLessThan(testCase,r.gradient_fd_relative_error,1e-6);
out=evalc("acoustic_fembem.check_energy_model_gate(61)"); j=jsondecode(out); verifyTrue(testCase,j.ok);
end


function testPolicyTrustGate(testCase)
r=acoustic_fembem.policy_trust_gate(67); verifyTrue(testCase,r.ok);
verifyLessThanOrEqual(testCase,r.trusted_kl,r.kl_limit); verifyLessThan(testCase,r.gradient_relative_error,1e-8);
out=evalc("acoustic_fembem.check_policy_trust_gate(67)"); j=jsondecode(out); verifyTrue(testCase,j.ok);
end


function testVariationalGeometryGate(testCase)
r=acoustic_fembem.variational_geometry_gate(71); verifyTrue(testCase,r.ok);
verifyGreaterThanOrEqual(testCase,r.kl,0); verifyLessThan(testCase,r.identity_error,1e-12);
verifyGreaterThan(testCase,r.bad_map_error,r.symplectic_error);
out=evalc("acoustic_fembem.check_variational_geometry_gate(71)"); j=jsondecode(out); verifyTrue(testCase,j.ok);
end


function testLatentDistributionGate(testCase)
r=acoustic_fembem.latent_distribution_gate(73); verifyTrue(testCase,r.ok);
verifyEqual(testCase,sum(r.posterior),1,"AbsTol",1e-12);
out=evalc("acoustic_fembem.check_latent_distribution_gate(73)"); j=jsondecode(out); verifyTrue(testCase,j.ok);
end


function testLinearInverseGate(testCase)
r=acoustic_fembem.linear_inverse_gate(79); verifyTrue(testCase,r.ok);
verifyLessThan(testCase,max(r.moore_penrose_errors),1e-10);
out=evalc("acoustic_fembem.check_linear_inverse_gate(79)"); j=jsondecode(out); verifyTrue(testCase,j.ok);
end


function testActorCriticGate(testCase)
r=acoustic_fembem.actor_critic_gate(83); verifyTrue(testCase,r.ok);
verifyLessThan(testCase,r.gradient_relative_error,1e-8);
verifyLessThan(testCase,r.baseline_variance_ratio,1);
out=evalc("acoustic_fembem.check_actor_critic_gate(83)"); j=jsondecode(out); verifyTrue(testCase,j.ok);
end


function testA3cAsyncGate(testCase)
r=acoustic_fembem.a3c_async_gate(233); verifyTrue(testCase,r.ok);
verifyTrue(testCase,r.checks.off_policy_negative_rejected);
verifyTrue(testCase,r.checks.stale_negative_rejected);
verifyGreaterThan(testCase,r.forward_objective_after,r.forward_objective_before);
verifyEqual(testCase,r.global_version_after-r.global_version_before,2);
r2=acoustic_fembem.a3c_async_gate(233);
verifyEqual(testCase,r.design_candidate,r2.design_candidate,"AbsTol",0);
out=evalc("acoustic_fembem.check_a3c_async_gate(233)"); j=jsondecode(out);
verifyTrue(testCase,j.ok); verifyEqual(testCase,j.tool,'acoustic_fembem_a3c_async_gate');
end


function testRandomMatrixStabilityGate(testCase)
r=acoustic_fembem.random_matrix_stability_gate(89); verifyTrue(testCase,r.ok);
verifyGreaterThan(testCase,r.exploding_max_singular,20);
verifyLessThan(testCase,r.vanishing_max_singular,.03);
out=evalc("acoustic_fembem.check_random_matrix_stability_gate(89)"); j=jsondecode(out); verifyTrue(testCase,j.ok);
end


function testBellmanControlGate(testCase)
r=acoustic_fembem.bellman_control_gate(97); verifyTrue(testCase,r.ok);
verifyLessThan(testCase,r.expectation_residual,1e-12);
verifyLessThan(testCase,r.value_iteration_final_residual,1e-4);
out=evalc("acoustic_fembem.check_bellman_control_gate(97)"); j=jsondecode(out); verifyTrue(testCase,j.ok);
end


function testSpectralLinearAlgebraGate(testCase)
r=acoustic_fembem.spectral_linear_algebra_gate(101); verifyTrue(testCase,r.ok);
verifyLessThan(testCase,r.svd_reconstruction_error,1e-12);
verifyLessThan(testCase,r.diagonalization_error,1e-12);
out=evalc("acoustic_fembem.check_spectral_linear_algebra_gate(101)"); j=jsondecode(out); verifyTrue(testCase,j.ok);
end


function testPcaCovarianceGate(testCase)
r=acoustic_fembem.pca_covariance_gate(103); verifyTrue(testCase,r.ok);
verifyLessThan(testCase,r.reconstruction_mse(3),1e-25);
out=evalc("acoustic_fembem.check_pca_covariance_gate(103)"); j=jsondecode(out); verifyTrue(testCase,j.ok);
end


function testQuantumLinearGate(testCase)
r=acoustic_fembem.quantum_linear_gate(107); verifyTrue(testCase,r.ok);
verifyLessThan(testCase,max(r.unitarity_errors),1e-12);
verifyEqual(testCase,r.deutsch_balanced_p1,1,"AbsTol",1e-12);
out=evalc("acoustic_fembem.check_quantum_linear_gate(107)"); j=jsondecode(out); verifyTrue(testCase,j.ok);
end


function testTransformerScalingGate(testCase)
r=acoustic_fembem.transformer_scaling_gate(109); verifyTrue(testCase,r.ok);
verifyGreaterThan(testCase,r.fit_r2,.99);
verifyFalse(testCase,r.extrapolation_verified);
out=evalc("acoustic_fembem.check_transformer_scaling_gate(109)"); j=jsondecode(out); verifyTrue(testCase,j.ok);
end


function testNlpTransformerGate(testCase)
r=acoustic_fembem.nlp_transformer_gate(113); verifyTrue(testCase,r.ok);
verifyLessThan(testCase,r.future_attention_leak,1e-12); verifyEqual(testCase,r.bleu_exact,1,"AbsTol",1e-12);
out=evalc("acoustic_fembem.check_nlp_transformer_gate(113)"); j=jsondecode(out); verifyTrue(testCase,j.ok);
end


function testQuantificationGate(testCase)
r=acoustic_fembem.quantification_gate(127); verifyTrue(testCase,r.ok);
verifyLessThan(testCase,r.reconstruction_error,1e-12);
out=evalc("acoustic_fembem.check_quantification_gate(127)"); j=jsondecode(out); verifyTrue(testCase,j.ok);
end


function testHypothesisTestGate(testCase)
r=acoustic_fembem.hypothesis_test_gate(131); verifyTrue(testCase,r.ok);
verifyGreaterThan(testCase,r.null_p,r.alpha); verifyLessThan(testCase,r.alternative_p,r.alpha);
out=evalc("acoustic_fembem.check_hypothesis_test_gate(131)"); j=jsondecode(out); verifyTrue(testCase,j.ok);
end


function testQuantificationIiGate(testCase)
r=acoustic_fembem.quantification_ii_gate(137); verifyTrue(testCase,r.ok);
verifyLessThan(testCase,r.lagrange_residual,1e-8);
out=evalc("acoustic_fembem.check_quantification_ii_gate(137)"); j=jsondecode(out); verifyTrue(testCase,j.ok);
end


function testFactorRotationGate(testCase)
r=acoustic_fembem.factor_rotation_gate(139); verifyTrue(testCase,r.ok);
verifyLessThan(testCase,r.promax_covariance_error,1e-12);
out=evalc("acoustic_fembem.check_factor_rotation_gate(139)"); j=jsondecode(out); verifyTrue(testCase,j.ok);
end


function testTextRetrievalGate(testCase)
r=acoustic_fembem.text_retrieval_gate(149); verifyTrue(testCase,r.ok);
verifyEqual(testCase,r.tfidf_scores(1),max(r.tfidf_scores));
out=evalc("acoustic_fembem.check_text_retrieval_gate(149)"); j=jsondecode(out); verifyTrue(testCase,j.ok);
end


function testProjectiveGeometryGate(testCase)
r=acoustic_fembem.projective_geometry_gate(151); verifyTrue(testCase,r.ok);
verifyLessThan(testCase,r.roundtrip_error,1e-12);
out=evalc("acoustic_fembem.check_projective_geometry_gate(151)"); j=jsondecode(out); verifyTrue(testCase,j.ok);
end


function testCnnArchitectureGate(testCase)
r=acoustic_fembem.cnn_architecture_gate(157); verifyTrue(testCase,r.ok);
verifyLessThan(testCase,r.separable_parameters,r.standard_parameters);
out=evalc("acoustic_fembem.check_cnn_architecture_gate(157)"); j=jsondecode(out); verifyTrue(testCase,j.ok);
end


function testRecurrentSequenceGate(testCase)
r=acoustic_fembem.recurrent_sequence_gate(163); verifyTrue(testCase,r.ok);
verifyGreaterThan(testCase,r.exploding_control_final,100);
out=evalc("acoustic_fembem.check_recurrent_sequence_gate(163)"); j=jsondecode(out); verifyTrue(testCase,j.ok);
end


function testLatentStructuralGate(testCase)
r=acoustic_fembem.latent_structural_gate(167); verifyTrue(testCase,r.ok);
verifyLessThan(testCase,r.posterior_group_mse,r.raw_group_mse);
out=evalc("acoustic_fembem.check_latent_structural_gate(167)"); j=jsondecode(out); verifyTrue(testCase,j.ok);
end


function testGraphicalModelGate(testCase)
r=acoustic_fembem.graphical_model_gate(173); verifyTrue(testCase,r.ok);
verifyLessThan(testCase,r.conditional_independence_error,1e-12);
out=evalc("acoustic_fembem.check_graphical_model_gate(173)"); j=jsondecode(out); verifyTrue(testCase,j.ok);
end


function testSdeItoGate(testCase)
r=acoustic_fembem.sde_ito_gate(179); verifyTrue(testCase,r.ok);
verifyLessThan(testCase,abs(r.pde_residual),1e-10);
out=evalc("acoustic_fembem.check_sde_ito_gate(179)"); j=jsondecode(out); verifyTrue(testCase,j.ok);
end


function testTimeSeriesGate(testCase)
r=acoustic_fembem.time_series_gate(181); verifyTrue(testCase,r.ok);
verifyLessThan(testCase,r.var_fit_error,.03);
out=evalc("acoustic_fembem.check_time_series_gate(181)"); j=jsondecode(out); verifyTrue(testCase,j.ok);
end


function testGradientBoostGate(testCase)
r=acoustic_fembem.gradient_boost_gate(191); verifyTrue(testCase,r.ok);
verifyLessThan(testCase,r.loss_history(end),r.loss_history(1));
out=evalc("acoustic_fembem.check_gradient_boost_gate(191)"); j=jsondecode(out); verifyTrue(testCase,j.ok);
end


function testNeuralOdeGate(testCase)
r=acoustic_fembem.neural_ode_gate(193); verifyTrue(testCase,r.ok);
verifyTrue(testCase,all(diff(r.errors)<0));
out=evalc("acoustic_fembem.check_neural_ode_gate(193)"); j=jsondecode(out); verifyTrue(testCase,j.ok);
end


function testTabularLookupGate(testCase)
r=acoustic_fembem.tabular_lookup_gate(197); verifyTrue(testCase,r.ok);
verifyEqual(testCase,r.exact_value,3.3);
out=evalc("acoustic_fembem.check_tabular_lookup_gate(197)"); j=jsondecode(out); verifyTrue(testCase,j.ok);
end


function testTabularTransformGate(testCase)
r=acoustic_fembem.tabular_transform_gate(199); verifyTrue(testCase,r.ok);
verifyEqual(testCase,r.sum,10);
out=evalc("acoustic_fembem.check_tabular_transform_gate(199)"); j=jsondecode(out); verifyTrue(testCase,j.ok);
end


function testBayesInferenceGate(testCase)
r=acoustic_fembem.bayes_inference_gate(211); verifyTrue(testCase,r.ok);
verifyLessThan(testCase,abs(r.mcmc_mean-r.posterior_mean),.015);
out=evalc("acoustic_fembem.check_bayes_inference_gate(211)"); j=jsondecode(out); verifyTrue(testCase,j.ok);
end


function testActivationProbabilityGate(testCase)
r=acoustic_fembem.activation_probability_gate(223); verifyTrue(testCase,r.ok);
verifyLessThan(testCase,r.softmax_jacobian_error,1e-9);
out=evalc("acoustic_fembem.check_activation_probability_gate(223)"); j=jsondecode(out); verifyTrue(testCase,j.ok);
end


function testTopologyInformationGate(testCase)
r=acoustic_fembem.topology_information_gate(227); verifyTrue(testCase,r.ok);
verifyTrue(testCase,all(diff(r.component_counts)<=0));
out=evalc("acoustic_fembem.check_topology_information_gate(227)"); j=jsondecode(out); verifyTrue(testCase,j.ok);
end


function testIntegrityEncodingGate(testCase)
r=acoustic_fembem.integrity_encoding_gate(229); verifyTrue(testCase,r.ok);
verifyNotEqual(testCase,r.hashes,r.tampered_hashes);
out=evalc("acoustic_fembem.check_integrity_encoding_gate(229)"); j=jsondecode(out); verifyTrue(testCase,j.ok);
end


function testKnowledgeIncludesPublicAcousticBlogLessons(testCase)
body = acoustic_fembem.fembem_knowledge("public_acoustic_blog_lessons");
verifyGreaterThan(testCase, strlength(body), 700);
verifySubstring(testCase, body, "Unbounded exterior radiation");
verifySubstring(testCase, body, "frequency domain");
verifySubstring(testCase, body, "high-order surface");
verifySubstring(testCase, body, "PML=false");
verifySubstring(testCase, body, "time-domain lane is explicitly CQ");
verifySubstring(testCase, body, "Acoustic-structure interaction");
verifySubstring(testCase, body, "two-way coupling");
verifySubstring(testCase, body, "Impedance lumping");
verifySubstring(testCase, body, "p=Z_s v");
verifySubstring(testCase, body, "p=Z Q");
verifySubstring(testCase, body, "local reaction and extended reaction");
verifySubstring(testCase, body, "high-order Zs");
verifySubstring(testCase, body, "Schroeder-frequency");
verifySubstring(testCase, body, "acoustic_method_selection_manifest_gate");
end


function testKnowledgeIncludesPublicAcousticNonboundary10(testCase)
body = acoustic_fembem.fembem_knowledge("public_acoustic_nonboundary_10");
verifyGreaterThan(testCase, strlength(body), 1200);
verifySubstring(testCase, body, "10 public acoustic non-boundary problems");
verifySubstring(testCase, body, "absorbing boundaries");
verifySubstring(testCase, body, "PML=false");
verifySubstring(testCase, body, "Acoustic trap");
verifySubstring(testCase, body, "Surface-acoustic-wave droplet streaming");
verifySubstring(testCase, body, "Thermoviscous acoustic radiation force");
verifySubstring(testCase, body, "Thermoacoustic engine");
verifySubstring(testCase, body, "Acoustic topology optimization");
verifySubstring(testCase, body, "Room response split");
verifySubstring(testCase, body, "Small-speaker room impulse response");
verifySubstring(testCase, body, "Ultrasonic pipe pulse-echo");
verifySubstring(testCase, body, "public_acoustic_nonboundary_problem_catalog");
verifySubstring(testCase, body, "acoustic_nonboundary_problem_catalog_manifest_gate");
end


function testKnowledgeIncludesGmshArtifactTopic(testCase)
% The MATLAB lane is gmsh-free: the topic documents the native GIF path and
% points gmsh acoustic movies to the radia-acoustic (GmshPostExport) side.
body = acoustic_fembem.fembem_knowledge("gmsh_artifact");
verifyGreaterThan(testCase, strlength(body), 500);
verifySubstring(testCase, body, "gmsh-free");
verifySubstring(testCase, body, "writeSoftSphereScatterGif");
verifySubstring(testCase, body, "drumScatterField");
verifySubstring(testCase, body, "does not require Gmsh");
verifySubstring(testCase, body, "GmshPostExport");
verifySubstring(testCase, body, "radia-acoustic");
verifySubstring(testCase, body, ".msh v4.1");
end


function testKnowledgeIncludesCatalog100Topic(testCase)
body = acoustic_fembem.fembem_knowledge("catalog_100");
verifyGreaterThan(testCase, strlength(body), 500);
verifySubstring(testCase, body, "100-case");
verifySubstring(testCase, body, "GYP-001..010");
verifySubstring(testCase, body, "GYP-091..100");
end


function testKnowledgeIncludesDrumTopic(testCase)
body = acoustic_fembem.fembem_knowledge("vibroacoustic_drum");
verifyGreaterThan(testCase, strlength(body), 500);
verifySubstring(testCase, body, "baffled circular membrane");
verifySubstring(testCase, body, "normal velocity");
verifySubstring(testCase, body, "NGSolve.BEM");
verifySubstring(testCase, body, "the drum structure is FEM");
verifySubstring(testCase, body, "air radiation");
verifySubstring(testCase, body, "acoustic BEM");
verifySubstring(testCase, body, "plotDrumStepTimeField");
verifySubstring(testCase, body, "writeDrumStepTimeGif");
verifySubstring(testCase, body, "drumHighOrderImpedanceScene");
verifySubstring(testCase, body, "drumFemBemCoupledDemo");
verifySubstring(testCase, body, "volFemBemIfftResponse");
verifySubstring(testCase, body, "volTdBemConvolutionQuadrature");
verifySubstring(testCase, body, "axis-equal");
verifySubstring(testCase, body, "not a hemisphere");
verifySubstring(testCase, body, "top membrane");
verifySubstring(testCase, body, "lower half-space is intentionally quiet");
verifySubstring(testCase, body, "rigid baffle");
verifySubstring(testCase, body, "high-order impedance boundary is mandatory");
verifySubstring(testCase, body, "do not use or name a Kelvin");
verifySubstring(testCase, body, "reduced FEM ODE");
verifySubstring(testCase, body, "ode45");
verifySubstring(testCase, body, "damping-ratio");
verifySubstring(testCase, body, "not a cavity");
verifySubstring(testCase, body, "decaying membrane/shell vibration");
verifySubstring(testCase, body, "same modeling split can be implemented in NGSolve");
verifySubstring(testCase, body, "retarded boundary");
verifySubstring(testCase, body, "must NOT split the observation field by source direction");
verifySubstring(testCase, body, "all evaluated at every");
verifySubstring(testCase, body, "direction-only painting is a");
verifySubstring(testCase, body, "not a cavity");
verifySubstring(testCase, body, "pressure DOF");
verifySubstring(testCase, body, "lower-half radiation");
verifySubstring(testCase, body, "3D axisymmetric");
verifySubstring(testCase, body, "r-z slice");
verifySubstring(testCase, body, "not yet a full 3D structural-FEM/acoustic-BEM drum mesh");
verifySubstring(testCase, body, "frequency-domain Helmholtz FEM/BEM");
verifySubstring(testCase, body, "inverse FFT");
verifySubstring(testCase, body, "not a periodic sine-wave animation");
verifySubstring(testCase, body, "parallel acoustic-volume teaching lane");
verifySubstring(testCase, body, "not the preferred");
verifySubstring(testCase, body, "BDF generating function");
verifySubstring(testCase, body, "Laplace-domain single-layer");
verifySubstring(testCase, body, "real Lubich CQ TD-BEM");
verifySubstring(testCase, body, "volFemBemCoupledConvolutionQuadrature");
verifySubstring(testCase, body, "H1/P1 interior wave FEM");
verifySubstring(testCase, body, "(1/2 Mb-K(s))*T");
verifySubstring(testCase, body, "-S(s)q + D(s)Tu");
verifySubstring(testCase, body, "Calderon/Johnson-Nedelec coupled CQ");
verifySubstring(testCase, body, "retarded double-layer K(s)");
verifySubstring(testCase, body, "SingleLayerTeaching");
verifySubstring(testCase, body, "GmshPostExport");
verifySubstring(testCase, body, "replace the interior acoustic volume FEM with structural membrane/shell FEM");
verifySubstring(testCase, body, "does not require Gmsh");
end


function testKnowledgeIncludesCurvedVolGeometryTopic(testCase)
body = acoustic_fembem.fembem_knowledge("curved_vol_geometry");
verifyGreaterThan(testCase, strlength(body), 500);
verifySubstring(testCase, body, "superparametric");
verifySubstring(testCase, body, "CurvedPanelQuadrature");
verifySubstring(testCase, body, "curvedSingleLayerDirichletSolve");
verifySubstring(testCase, body, "curvedelements");
end


function testKnowledgeIncludesMatlabExecutionPolicy(testCase)
body = acoustic_fembem.fembem_knowledge("matlab_execution_policy");
verifyGreaterThan(testCase, strlength(body), 300);
verifySubstring(testCase, body, ".m functions/scripts");
verifySubstring(testCase, body, "MCP tools");
verifySubstring(testCase, body, "JSON manifests");
end


function testKnowledgeIncludesMathWorksAgenticToolkitPolicy(testCase)
body = acoustic_fembem.fembem_knowledge("mathworks_agentic_toolkit");
verifyGreaterThan(testCase, strlength(body), 500);
verifySubstring(testCase, body, "official MathWorks MATLAB MCP Server");
verifySubstring(testCase, body, "MATLAB Agentic Toolkit");
verifySubstring(testCase, body, "runtime");
verifySubstring(testCase, body, "skills");
verifySubstring(testCase, body, "existing-session");
verifySubstring(testCase, body, "acoustic_fembem extension");
verifyEqual(testCase, acoustic_fembem.fembem_knowledge("matlab_mcp_server"), body);
end


function testVolMeshSummaryWrapper(testCase)
out = evalc("acoustic_fembem.check_vol_mesh_summary(""unit_sphere_coarse.vol"")");
decoded = jsondecode(out);
verifyTrue(testCase, decoded.ok);
verifyEqual(testCase, string(decoded.tool), "acoustic_fembem_vol_mesh_summary");
verifyEqual(testCase, string(decoded.recommended_gui_viewer), "Netgen/native .vol viewer");
verifyTrue(testCase, contains(string(decoded.recommended_windows_double_click_handler), "Ng_LoadMesh"));
verifyGreaterThan(testCase, decoded.points, 0);
verifyGreaterThan(testCase, decoded.triangles, 0);
verifyGreaterThan(testCase, decoded.tets, 0);
end


function testRepositoryHealthWrapper(testCase)
out = evalc("acoustic_fembem.check_repository_health()");
decoded = jsondecode(out);
verifyTrue(testCase, decoded.ok);
verifyEqual(testCase, string(decoded.tool), "acoustic_fembem_repository_health");
verifyEqual(testCase, string(decoded.repository_name), "matlab-acoustic-fembem");
verifyEqual(testCase, decoded.num_validation_cases, 100);
verifyEqual(testCase, decoded.num_verified_cases, 100);
verifyGreaterThanOrEqual(testCase, decoded.num_vol_fixtures, 10);
end


function testResultManifestGateWrapper(testCase)
artifact = completeArtifact();
manifestPath = fullfile(tempdir, "acoustic_fembem_result_manifest_gate_test.json");
fid = fopen(manifestPath, "w");
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, "%s", jsonencode(artifact));
clear cleanup

out = evalc("acoustic_fembem.check_result_manifest_file(manifestPath, true, true, ""matlab_fem_bem_result_table_v1"")");
decoded = jsondecode(out);
verifyTrue(testCase, decoded.ok);
verifyEqual(testCase, string(decoded.tool), "acoustic_fembem_check_result_manifest_file");
end


function testAnalyticAcousticGateWrapper(testCase)
out = evalc("acoustic_fembem.check_fembem_acoustic_gate(""soft"", 2.0, 7, -1)");
decoded = jsondecode(out);
verifyTrue(testCase, decoded.ok);
verifyEqual(testCase, string(decoded.tool), "acoustic_fembem_acoustic_gate");
verifyEqual(testCase, string(decoded.kind), "soft");
end


function testCrossvalGateWrapper(testCase)
out = evalc("acoustic_fembem.check_fembem_crossval_gate(""galerkin_ngsolve"", ""unit_sphere_coarse.vol"", -1, false)");
decoded = jsondecode(out);
verifyTrue(testCase, decoded.ok);
verifyEqual(testCase, string(decoded.tool), "acoustic_fembem_crossval_gate");
verifyEqual(testCase, string(decoded.input_format), "netgen_vol_tri_tet");
end


function artifact = completeArtifact()
artifact = struct();
artifact.schema = "matlab-acoustic-fembem.crossval.v1";
artifact.pass = true;
artifact.created_at_utc = "2026-07-04T00:00:00Z";
artifact.versions = struct("matlab", version, "radia_mcp", "test");
artifact.execution = struct( ...
    "run_date_utc", "2026-07-04T00:00:02Z", ...
    "execution_session_id", "MATLAB_TEST");
artifact.expected_created_at_utc = "2026-07-04T00:00:00Z";
artifact.expected_run_date_utc = "2026-07-04T00:00:02Z";
artifact.expected_execution_session_id = "MATLAB_TEST";
artifact.result_output_schema_id = "matlab_fem_bem_result_table_v1";
artifact.result_output_columns = ["alpha", "trace_residual_norm"];
artifact.result_output_units = struct("alpha", "1", "trace_residual_norm", "1");
artifact.timing_breakdown_s = struct("solve", 0.1, "postprocess", 0.02);
artifact.physics_convention_schema_id = "matlab_first_order_fem_bem_coupling_convention_v1";
artifact.postprocess_row_convention_schema_id = "matlab_fem_bem_postprocess_rows_v1";
end
