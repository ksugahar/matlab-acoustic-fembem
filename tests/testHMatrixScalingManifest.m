function tests = testHMatrixScalingManifest
tests = functiontests(localfunctions);
end


function setupOnce(~)
repoRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(repoRoot);
addpath(genpath(fullfile(repoRoot, "matlab_api")));
end


function testStudyShowsBoundedRankAndLinearStorage(testCase)
study = hmatrixScalingStudy([60 120 240]);
report = hmatrixScalingManifest(study.rows);
verifyTrue(testCase, report.ok);
verifyEqual(testCase, report.maxRanks, [5 5 5]);
verifyEqual(testCase, report.storageGrowthExponents, [1 1], "AbsTol", 1e-12);
verifyLessThan(testCase, report.maxMatvecRelativeError, 1e-10);
end


function testManifestRejectsRankAndStorageRegression(testCase)
study = hmatrixScalingStudy([30 60 120]);
rows = study.rows;
rows(3).maxRank = 40;
rows(3).storedEntries = 10000;
rows(3).compressionRatio = rows(3).storedEntries/rows(3).denseEntries;
report = hmatrixScalingManifest(rows);
verifyFalse(testCase, report.ok);
verifyFalse(testCase, report.checks.rankIsBounded);
verifyFalse(testCase, report.checks.storedEntryGrowthIsSubquadratic);
end


function testMcpWrapperRunsPositiveAndNegativeThresholds(testCase)
out = evalc("acoustic_fembem.check_hmatrix_scaling(""[60,120,240]"", 1e-8, 1e-8, 20, 1.25)");
decoded = jsondecode(out);
verifyTrue(testCase, decoded.ok);
verifyEqual(testCase, string(decoded.tool), "acoustic_fembem_hmatrix_scaling");

bad = evalc("acoustic_fembem.check_hmatrix_scaling(""[60,120,240]"", 1e-8, 1e-8, 4, 1.25)");
rejected = jsondecode(bad);
verifyFalse(testCase, rejected.ok);
verifyFalse(testCase, rejected.result.checks.rankIsBounded);
end
