@echo off
setlocal EnableExtensions DisableDelayedExpansion

if /i "%~1"=="--help" goto :help
if /i "%~1"=="/?" goto :help
if not "%~2"=="" goto :invalid_arguments
set "EC_BUILD_RUNTIME=win-x64"
if "%~1"=="" goto :start
if /i "%~1"=="win-x64" goto :start
if /i not "%~1"=="win-arm64" goto :invalid_arguments
set "EC_BUILD_RUNTIME=win-arm64"

:start
pushd "%~dp0"
if errorlevel 1 exit /b 1
set "EC_BUILD_ROOT=%CD%"

where.exe dotnet >nul 2>&1
if errorlevel 1 (
    echo ERROR: Install the .NET 8 SDK and make dotnet available on PATH.
    goto :failed
)
where.exe powershell.exe >nul 2>&1
if errorlevel 1 (
    echo ERROR: Windows PowerShell is required to create the release archives.
    goto :failed
)

set "EC_BUILD_STAMP="
for /f %%I in ('powershell.exe -NoProfile -NonInteractive -Command "Get-Date -Format yyyyMMdd-HHmmss-fff"') do set "EC_BUILD_STAMP=%%I"
if not defined EC_BUILD_STAMP goto :failed
set "EC_BUILD_OUTPUT=%EC_BUILD_ROOT%\artifacts\release\%EC_BUILD_STAMP%-%EC_BUILD_RUNTIME%"
if exist "%EC_BUILD_OUTPUT%" (
    echo ERROR: The release output directory already exists. Run the script again.
    goto :failed
)
mkdir "%EC_BUILD_OUTPUT%"
if errorlevel 1 goto :failed

echo [1/7] Building the solution in Release...
dotnet build EllipticCurves.sln -c Release -p:GeneratePackageOnBuild=false
if errorlevel 1 goto :failed

echo [2/7] Running arithmetic, Console and portable Explorer tests...
dotnet test tests\EllipticCurves.Tests.csproj -c Release --no-build --no-restore
if errorlevel 1 goto :failed

echo [3/7] Building and running the WPF presentation check...
dotnet build tests\PresentationHost\PresentationHost.csproj -c Release -p:GeneratePackageOnBuild=false
if errorlevel 1 goto :failed
dotnet run --project tests\PresentationHost\PresentationHost.csproj -c Release --no-build --no-restore
if errorlevel 1 goto :failed

echo [4/7] Packing the library...
dotnet pack sources\EllipticCurves.csproj -c Release --no-build --no-restore -p:GeneratePackageOnBuild=false -o "%EC_BUILD_OUTPUT%\nuget"
if errorlevel 1 goto :failed

echo [5/7] Publishing Explorer with the .NET runtime for %EC_BUILD_RUNTIME%...
dotnet publish explorer\EllipticCurves.Explorer.csproj -c Release -r %EC_BUILD_RUNTIME% --self-contained true -p:GeneratePackageOnBuild=false -o "%EC_BUILD_OUTPUT%\explorer-%EC_BUILD_RUNTIME%"
if errorlevel 1 goto :failed

echo [6/7] Publishing Console with the .NET runtime for %EC_BUILD_RUNTIME%...
dotnet publish console\EllipticCurves.Console.csproj -c Release -r %EC_BUILD_RUNTIME% --self-contained true -p:GeneratePackageOnBuild=false -o "%EC_BUILD_OUTPUT%\console-%EC_BUILD_RUNTIME%"
if errorlevel 1 goto :failed

echo [7/7] Adding licenses and creating ZIP archives...
powershell.exe -NoProfile -NonInteractive -Command "$ErrorActionPreference = 'Stop'; $version = & dotnet msbuild (Join-Path $env:EC_BUILD_ROOT 'sources/EllipticCurves.csproj') -nologo -p:Configuration=Release -getProperty:PackageVersion; if ($LASTEXITCODE -ne 0 -or [string]::IsNullOrWhiteSpace($version)) { throw 'Cannot read the library package version.' }; $license = Join-Path $env:EC_BUILD_ROOT 'LICENSE'; $explorer = Join-Path $env:EC_BUILD_OUTPUT ('explorer-' + $env:EC_BUILD_RUNTIME); $console = Join-Path $env:EC_BUILD_OUTPUT ('console-' + $env:EC_BUILD_RUNTIME); foreach ($folder in @($explorer, $console)) { Copy-Item -LiteralPath $license -Destination (Join-Path $folder 'EllipticCurves.LICENSE.txt') }; Compress-Archive -Path (Join-Path $explorer '*') -DestinationPath (Join-Path $env:EC_BUILD_OUTPUT ('EllipticCurves.Explorer.' + $version + '-' + $env:EC_BUILD_RUNTIME + '.zip')); Compress-Archive -Path (Join-Path $console '*') -DestinationPath (Join-Path $env:EC_BUILD_OUTPUT ('EllipticCurves.Console.' + $version + '-' + $env:EC_BUILD_RUNTIME + '.zip'))"
if errorlevel 1 goto :failed

echo.
echo Release artifacts are ready:
echo "%EC_BUILD_OUTPUT%"
popd
endlocal
exit /b 0

:failed
set "EC_BUILD_EXIT=%errorlevel%"
if "%EC_BUILD_EXIT%"=="0" set "EC_BUILD_EXIT=1"
echo.
echo ERROR: Release preparation failed. See the output above.
popd
endlocal & exit /b %EC_BUILD_EXIT%

:invalid_arguments
echo ERROR: Expected win-x64 or win-arm64.
echo Usage: build.bat [win-x64 ^| win-arm64]
exit /b 2

:help
echo Usage: build.bat [win-x64 ^| win-arm64]
echo Default: win-x64. Requires the .NET 8 SDK and Windows PowerShell.
echo Builds Release, runs tests, and creates NuGet, Explorer and Console artifacts.
echo Explorer and Console include .NET and run without a separate runtime installation.
echo Output: artifacts\release\TIMESTAMP-RUNTIME
echo Artifacts are created locally; nothing is uploaded.
exit /b 0
