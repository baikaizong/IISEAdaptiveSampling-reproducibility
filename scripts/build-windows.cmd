@echo off
setlocal EnableExtensions EnableDelayedExpansion

for %%I in ("%~dp0..") do set "REPOSITORY_ROOT=%%~fI"
if not defined CMAKE_EXE set "CMAKE_EXE=cmake"

"%CMAKE_EXE%" --version >nul 2>&1
if errorlevel 1 (
  echo CMake is not available. Add it to PATH or set CMAKE_EXE.
  exit /b 2
)

pushd "%REPOSITORY_ROOT%"
"%CMAKE_EXE%" --build --preset windows-debug --parallel
set "EXIT_CODE=!ERRORLEVEL!"
popd
exit /b !EXIT_CODE!
