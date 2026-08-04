@echo off
setlocal EnableExtensions EnableDelayedExpansion

for %%I in ("%~dp0..") do set "REPOSITORY_ROOT=%%~fI"

if "%~1"=="" goto usage
set "APP=%~1"
set "EXPERIMENT=%~2"
set "NPROCS=%~3"
set "RESULTS_ROOT=%~4"
set "CACHE_ROOT=%~5"
set "SEED=%~6"

if "%EXPERIMENT%"=="" set "EXPERIMENT=default"
if "%NPROCS%"=="" set "NPROCS=4"
if "%RESULTS_ROOT%"=="" set "RESULTS_ROOT=results"
if "%CACHE_ROOT%"=="" set "CACHE_ROOT=cache"
if "%SEED%"=="" set "SEED=20260804"

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

if "%EXPERIMENT%"=="default" (
  if /I "%APP%"=="mmkstd" set "EXPERIMENT=circle-comparison"
  if /I "%APP%"=="cds" set "EXPERIMENT=kernel-circle"
  if /I "%APP%"=="nas" set "EXPERIMENT=kernel-circle"
  if /I "%APP%"=="pos" set "EXPERIMENT=kernel-circle"
  if /I "%APP%"=="sasam" set "EXPERIMENT=kernel-circle"
  if /I "%APP%"=="topr" set "EXPERIMENT=kernel-circle"
  if /I "%APP%"=="tsbss" set "EXPERIMENT=bspline"
  if /I "%APP%"=="tsspr" set "EXPERIMENT=kernel-circle"
)

set "EXE=build\windows-debug\bin\%APP%_evaluation.exe"
pushd "%REPOSITORY_ROOT%"
if not exist "%EXE%" (
  echo Executable not found: !EXE!
  popd
  exit /b 2
)
mpiexec -genv PATH "!PATH!" -n %NPROCS% "%EXE%" "%EXPERIMENT%" "%RESULTS_ROOT%" "%CACHE_ROOT%" "%SEED%"
set "EXIT_CODE=!ERRORLEVEL!"
popd
exit /b !EXIT_CODE!

:usage
echo Usage: scripts\run-windows.cmd APP [EXPERIMENT] [MPI_PROCESSES] [RESULTS_ROOT] [CACHE_ROOT] [SEED]
echo APP: mmkstd, cds, nas, pos, sasam, topr, tsbss, or tsspr
exit /b 2
