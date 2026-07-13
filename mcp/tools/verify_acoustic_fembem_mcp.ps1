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

$batchCode = "repoRoot='$($repoRoot.Path.Replace("'", "''"))'; addpath(repoRoot); addpath(genpath(fullfile(repoRoot,'matlab_api'))); addpath(fullfile(repoRoot,'examples')); addpath(fullfile(repoRoot,'validation')); results=runtests(fullfile(repoRoot,'tests','testMcpAcousticFembemTools.m')); assertSuccess(results);"
if ($SkipSharedEngine) {
    Write-Warning "Shared engine checks skipped (-SkipSharedEngine). Using matlab -batch for pure MATLAB validation."
    matlab -batch $batchCode
    exit $LASTEXITCODE
}

$sharedEngineCode = @'
import sys
import matlab.engine

names = matlab.engine.find_matlab()
print(f"[acoustic-fembem-mcp] shared engines: {names}")
if not names:
    raise SystemExit(2)
engine = matlab.engine.connect_matlab(names[0])
try:
    output = engine.evalc(sys.argv[1])
    if output:
        print(output, end="")
except BaseException as exc:
    print(f"{type(exc).__name__}: {exc}", file=sys.stderr)
    raise SystemExit(1)
'@
$sharedEngineCode | python - $batchCode
exit $LASTEXITCODE
