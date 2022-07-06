@echo off
rem tch.bat

rem This is a DOS batch file to run H2OI95 test cases.
rem It compares present results with previous results.
rem It uses two other batch files, tchs.bat and tfcm1.bat.
rem It does not use tfcm1.bat directly, but through
rem tchs.bat.

rem Last revised 05/18/2020 by TJW

setlocal enabledelayedexpansion enableextensions

echo Running tch.bat-

rem Note: in setting codeDIR, it may be necessary to change
rem the drive letter.

set codeDIR=C:\H2OI95_v1.1

set codeEXE=%codeDIR%\Bin\H2OI95.exe
if exist %codeEXE% goto OK0
  echo     Error -(%0): Code executable "%codeEXE%" is not present.
  goto QUIT
:OK0

set bfserr=FALSE

set batfs1=%codeDIR%\Testing\tchs.bat
if exist %batfs1% goto OK1
  echo   Error - (%0): Batch file %batfs1% is missing.
  set bfserr=TRUE
:OK1

set batfs2=%codeDIR%\Testing\tfcm1.bat
if exist %batfs2% goto OK2
  echo   Error - (%0): Batch file %batfs2% is missing.
  set bfserr=TRUE
:OK2

rem If either subsidiary batch file is missing, quit.

if %bfserr%==TRUE goto QUIT

rem There should be one argument, the desired test
rem subdirectory.

if %1.==. goto SKIP1
if NOT %2.==. goto SKIP1

set testdir=%1
echo Test subdirectory: %1

set baseDIR=%codeDIR%\Testing\%testdir%
if NOT exist %baseDIR% goto SKIP2

echo Base directory is: %baseDIR%

cd %baseDIR%

for /f "eol=: delims=" %%F in ('dir /b /ad *') do (
  echo\
  echo Found folder %%F
  CALL %batfs1% %%F
  echo\
)

goto QUIT

:SKIP0
echo     Error -(%0): Code executable "%codeEXE%" is not present.
goto QUIT

:SKIP1
echo Usage: %0 testsubdirectory
goto QUIT

:SKIP2
echo   Error - (%0): Base directory %baseDIR% is missing.
goto QUIT

:QUIT
cd %codeDIR%\Testing
endlocal
echo All done
