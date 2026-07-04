function check_repository_health()
%CHECK_REPOSITORY_HEALTH Print a compact JSON verdict for repository health.
%
% MCP-facing entry point for a lightweight pre-push smoke check.

report = acoustic_fembem.repository_health();

summary = struct();
summary.tool = report.tool;
summary.status = report.status;
summary.ok = report.pass;
summary.repository_name = report.repository_name;
summary.repository_url = report.repository_url;
summary.num_validation_cases = report.num_validation_cases;
summary.num_verified_cases = report.num_verified_cases;
summary.num_vol_fixtures = report.num_vol_fixtures;
summary.failed_checks = report.failed_checks;

disp(jsonencode(summary));

if summary.ok
    return
end

error("acoustic_fembem:RepositoryHealthNeedsAttention", ...
    "Repository health needs attention. Failed checks: %s", ...
    strjoin(summary.failed_checks, ", "));
end
