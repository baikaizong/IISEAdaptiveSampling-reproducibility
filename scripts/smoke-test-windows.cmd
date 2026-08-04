@echo off
setlocal EnableExtensions EnableDelayedExpansion

for %%I in ("%~dp0..") do set "REPOSITORY_ROOT=%%~fI"
if not defined GUROBI_ROOT if defined GUROBI_HOME set "GUROBI_ROOT=%GUROBI_HOME%"
if not defined IMSL_ROOT (
  echo IMSL_ROOT is not set.
  exit /b 2
)
if not defined GUROBI_ROOT (
  echo GUROBI_ROOT is not set and GUROBI_HOME is unavailable.
  exit /b 2
)
if not exist "%IMSL_ROOT%\lib\imsl_dll.lib" (
  echo IMSL runtime library not found under: !IMSL_ROOT!
  exit /b 2
)
if not exist "%GUROBI_ROOT%\bin\gurobi130.dll" (
  echo Gurobi runtime not found under: !GUROBI_ROOT!
  exit /b 2
)
where mpiexec >nul 2>&1
if errorlevel 1 (
  echo mpiexec is not available. Initialize the Intel MPI environment first.
  exit /b 2
)
set "PATH=%IMSL_ROOT%\lib;%GUROBI_ROOT%\bin;%PATH%"

pushd "%REPOSITORY_ROOT%"
for %%A in (mmkstd cds nas pos sasam topr tsbss tsspr) do (
  echo Smoke testing %%A_evaluation.exe
  mpiexec -genv PATH "!PATH!" -n 2 "build\windows-debug\bin\%%A_evaluation.exe" invalid ^
    "build\smoke\results" "build\smoke\cache" 12345
  set "EXIT_CODE=!ERRORLEVEL!"
  if not "!EXIT_CODE!"=="0" (
    popd
    exit /b !EXIT_CODE!
  )
)
popd
echo All executable startup checks passed.
