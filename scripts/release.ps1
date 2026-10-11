<#
.SYNOPSIS
    Tags the current commit with the package version and pushes the tag, which publishes
    the package to NuGet.

.DESCRIPTION
    Reads <Version> from CSharpNumerics.csproj and creates the matching v<version> tag, so the
    tag and the package version can never disagree. Pushing a v* tag starts publish.yml, which
    runs dotnet pack and dotnet nuget push against nuget.org.

    That is a public, irreversible action: a version pushed to NuGet cannot be replaced, only
    deprecated or unlisted. So the script refuses to run unless every precondition holds, and
    asks for confirmation before the push.

    Checks, in order:
      1. The csproj declares a version, and it parses as a release version.
      2. The working tree is clean, so the tag names exactly what was tested.
      3. HEAD is on the release branch (master by default).
      4. HEAD exists on the remote, so the commit being published is one others can see.
      5. The tag does not already exist locally or on the remote.
      6. Optionally, the test suite passes (-Test).

.PARAMETER Branch
    Branch releases must come from. Defaults to master.

.PARAMETER Test
    Run the test suite in Release before tagging. Takes a couple of minutes.

.PARAMETER DryRun
    Run every check and print what would happen, without creating or pushing anything.

.PARAMETER Force
    Skip the confirmation prompt. The checks still run. Intended for CI, not for convenience.

.EXAMPLE
    ./scripts/release.ps1 -DryRun
    Verify that everything is in order for a release without touching anything.

.EXAMPLE
    ./scripts/release.ps1 -Test
    Run the suite, then tag and push after confirming.

.NOTES
    To release a new version: set <Version> in CSharpNumerics.csproj, commit, push, then run
    this. The script deliberately does not edit the version itself — bumping it is a decision
    that belongs in a reviewable commit, not in a release command.
#>
[CmdletBinding()]
param(
    [string]$Branch = 'master',
    [switch]$Test,
    [switch]$DryRun,
    [switch]$Force
)

$ErrorActionPreference = 'Stop'

$repoRoot = Split-Path -Parent $PSScriptRoot
$csproj = Join-Path $repoRoot 'Numerics/Numerics/CSharpNumerics.csproj'

function Fail($message) {
    Write-Host "FAILED: $message" -ForegroundColor Red
    exit 1
}

function Step($message) {
    Write-Host "  $message" -ForegroundColor DarkGray
}

Write-Host ''
Write-Host 'Release checks' -ForegroundColor Cyan

# ── 1. Version ────────────────────────────────────────────────────────────────
if (-not (Test-Path $csproj)) { Fail "Project file not found at $csproj" }

$versionMatch = Select-String -Path $csproj -Pattern '<Version>(.+?)</Version>' | Select-Object -First 1
if (-not $versionMatch) { Fail 'No <Version> element in CSharpNumerics.csproj' }

$version = $versionMatch.Matches[0].Groups[1].Value.Trim()
if ($version -notmatch '^\d+\.\d+\.\d+$') {
    Fail "Version '$version' is not a three-part release version. Pre-release suffixes need a deliberate decision about publishing them."
}

$tag = "v$version"
Step "version $version from csproj, tag $tag"

# ── 2. Clean working tree ─────────────────────────────────────────────────────
$dirty = git status --porcelain
if ($dirty) {
    Write-Host $dirty
    Fail 'Working tree has uncommitted changes. The tag must name exactly what was tested.'
}
Step 'working tree clean'

# ── 3. Right branch ───────────────────────────────────────────────────────────
$current = (git rev-parse --abbrev-ref HEAD).Trim()
if ($current -ne $Branch) {
    Fail "On branch '$current', but releases come from '$Branch'. Pass -Branch to override."
}
Step "on $Branch"

# ── 4. HEAD is on the remote ──────────────────────────────────────────────────
git fetch --quiet origin $Branch 2>$null
if ($LASTEXITCODE -ne 0) { Fail "Could not reach origin to verify $Branch." }

$head = (git rev-parse HEAD).Trim()
$remoteHasHead = git branch --remotes --contains $head 2>$null | Select-String "origin/$Branch"
if (-not $remoteHasHead) {
    Fail "HEAD ($($head.Substring(0,8))) is not on origin/$Branch. Push the branch first, so the published commit is one others can see."
}
Step 'HEAD present on origin'

# ── 5. Tag is free ────────────────────────────────────────────────────────────
$localTag = git tag --list $tag
if ($localTag) { Fail "Tag $tag already exists locally. Bump <Version> for a new release." }

$remoteTag = git ls-remote --tags origin "refs/tags/$tag" 2>$null
if ($remoteTag) { Fail "Tag $tag already exists on origin. Version $version has been released." }
Step "$tag is free"

# ── 6. Tests ──────────────────────────────────────────────────────────────────
if ($Test) {
    Write-Host '  running tests (Release)...' -ForegroundColor DarkGray
    $testProject = Join-Path $repoRoot 'Numerics/NumericTest/NumericTest.csproj'
    dotnet test $testProject -c Release --nologo
    if ($LASTEXITCODE -ne 0) { Fail 'Tests failed. Not tagging.' }
    Step 'tests passed'
}
else {
    Step 'tests skipped (pass -Test to run them)'
}

# ── Summary and confirmation ──────────────────────────────────────────────────
$subject = (git log -1 --pretty=%s).Trim()

Write-Host ''
Write-Host 'Ready to release' -ForegroundColor Cyan
Write-Host "  tag      $tag"
Write-Host "  commit   $($head.Substring(0,8))  $subject"
Write-Host "  effect   pushing $tag starts publish.yml, which pushes CSharpNumerics $version to nuget.org"
Write-Host ''

if ($DryRun) {
    Write-Host 'Dry run: nothing was created or pushed.' -ForegroundColor Yellow
    exit 0
}

if (-not $Force) {
    Write-Host 'This publishes publicly and cannot be undone.' -ForegroundColor Yellow
    $answer = Read-Host "Type $version to confirm"
    if ($answer -ne $version) {
        Write-Host 'Cancelled.' -ForegroundColor Yellow
        exit 1
    }
}

# ── Tag and push ──────────────────────────────────────────────────────────────
git tag -a $tag -m "Release $version"
if ($LASTEXITCODE -ne 0) { Fail 'Could not create the tag.' }
Write-Host "Created $tag" -ForegroundColor Green

git push origin $tag
if ($LASTEXITCODE -ne 0) {
    # Leave no half-finished state: a local tag that was never pushed is a trap for the
    # next run, which would refuse because the tag "already exists".
    git tag -d $tag | Out-Null
    Fail "Could not push the tag; removed it locally so this can be retried."
}

Write-Host "Pushed $tag" -ForegroundColor Green
Write-Host ''
Write-Host 'publish.yml is now running:' -ForegroundColor Cyan
Write-Host '  https://github.com/backlundtransform/CSharpNumerics/actions/workflows/publish.yml'
Write-Host ''
