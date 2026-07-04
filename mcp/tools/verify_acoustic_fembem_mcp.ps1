param(
    [switch]$SkipSharedEngine
)

$ErrorActionPreference = "Stop"

$repoRoot = Resolve-Path (Join-Path $PSScriptRoot "..\..")
$extensionFile = Join-Path $repoRoot "mcp\extensions\acoustic-fembem-tools.json"
$testFile = Join-Path $repoRoot "tests\testMcpAcousticFembemTools.m"

Write-Host "[acoustic-fembem-mcp] repo: $repoRoot"

if (-not (Test-Path -LiteralPath $extensionFile)) {
    throw "Missing extension file: $extensionFile"
}

if (-not (Test-Path -LiteralPath $testFile)) {
    throw "Missing test file: $testFile"
}

Get-Content -LiteralPath $extensionFile -Raw | ConvertFrom-Json | Out-Null

if (-not $SkipSharedEngine) {
    $findCode = @'
import matlab.engine
names = matlab.engine.find_matlab()
print(names)
raise SystemExit(0 if names else 2)
'@
    $findCode | python -
}

if ($SkipSharedEngine) {
    Write-Warning "Shared engine checks skipped (-SkipSharedEngine). Falling back to matlab -batch."
} else {
    Write-Warning "Existing MATLAB Engine discovery completed. Falling back to matlab -batch for test execution."
}

$batchCode = "repoRoot='$($repoRoot.Path.Replace("'", "''"))'; addpath(repoRoot); addpath(genpath(fullfile(repoRoot,'matlab_api'))); addpath(fullfile(repoRoot,'examples')); addpath(fullfile(repoRoot,'validation')); results=runtests(fullfile(repoRoot,'tests','testMcpAcousticFembemTools.m')); assertSuccess(results);"
matlab -batch $batchCode
exit $LASTEXITCODE
