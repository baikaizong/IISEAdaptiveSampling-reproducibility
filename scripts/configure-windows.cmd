@echo off
setlocal EnableExtensions EnableDelayedExpansion

for %%I in ("%~dp0..") do set "REPOSITORY_ROOT=%%~fI"

if not "%~1"=="" set "IMSL_ROOT=%~1"
if not "%~2"=="" set "GUROBI_ROOT=%~2"
if not defined GUROBI_ROOT if defined GUROBI_HOME set "GUROBI_ROOT=%GUROBI_HOME%"
if not defined CMAKE_EXE set "CMAKE_EXE=cmake"

if not defined IMSL_ROOT (
  echo IMSL_ROOT is not set.
  echo Set it to the architecture-specific IMSL installation directory.
  exit /b 2
)
if not defined GUROBI_ROOT (
  echo GUROBI_ROOT is not set and GUROBI_HOME is unavailable.
  echo Set either variable to the Gurobi installation directory.
  exit /b 2
)
where ifx >nul 2>&1
if errorlevel 1 (
  echo ifx is not available. Initialize the Intel oneAPI x64 environment first.
  exit /b 2
)
where icx >nul 2>&1
if errorlevel 1 (
  echo icx is not available. Initialize the Intel oneAPI x64 environment first.
  exit /b 2
)
where mpiexec >nul 2>&1
if errorlevel 1 (
  echo mpiexec is not available. Initialize the Intel MPI environment first.
  exit /b 2
)
where ninja >nul 2>&1
if errorlevel 1 (
  echo Ninja is not available on PATH.
  exit /b 2
)
"%CMAKE_EXE%" --version >nul 2>&1
if errorlevel 1 (
  echo CMake is not available. Add it to PATH or set CMAKE_EXE.
  exit /b 2
)
if not exist "%IMSL_ROOT%\include\dll\RNSET_INT.mod" (
  echo IMSL x64 DLL interface not found under: !IMSL_ROOT!
  exit /b 2
)
if not exist "%GUROBI_ROOT%\include\gurobi_c.h" (
  echo Gurobi headers not found under: !GUROBI_ROOT!
  exit /b 2
)

pushd "%REPOSITORY_ROOT%"
"%CMAKE_EXE%" --fresh --preset windows-debug -DIMSL_ROOT="%IMSL_ROOT%" -DGUROBI_ROOT="%GUROBI_ROOT%"
set "EXIT_CODE=!ERRORLEVEL!"
popd
exit /b !EXIT_CODE!
