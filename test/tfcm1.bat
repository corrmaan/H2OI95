@echo off
rem tfcm1.bat

rem Compare two text files, omitting the first line of each, which is presumed
rem to contain version numbering. A version of the sed editor is required to
rem make a copy of a file, omitting the first line. This batch file is intended
rem to support tchs.bat, which in turn supports tch.bat.

rem Usage: tfcm1 file1 file2

rem Last revised 04/18/2020 by TJW

rem echo   Executing tfcm1 %1 %2

if %1.==. goto SKIP1

if %2.==. goto SKIP2

if NOT %3.==. goto SKIP3

if NOT exist %1 goto NO1

if NOT exist %2 goto NO2

if exist acopy1.tmp del acopy1.tmp
if exist acopy2.tmp del acopy2.tmp

rem Make a copy of each file, deleting the first line of each.

sed '1d' %1 > acopy1.tmp
sed '1d' %2 > acopy2.tmp

set errstr=FALSE

fc acopy1.tmp acopy2.tmp > nul

if %errorlevel% == 0 goto OKAY
  echo   The files %1 and %2 are different

rem The value of errstr set here is effectively "returned"
rem to tchs.bat. The variable errorlevel apparently only
rem has local scope here. The value determined here is not
rem "returned" to tchs.bat.

  set errstr=TRUE

  goto ENDOKAY
:OKAY
  echo   The files %1 and %2 are identical
:ENDOKAY

del acopy1.tmp
del acopy2.tmp

goto QUIT

:SKIP1
echo     Error - (%0): First filename is missing. Usage: tfcm1 file1 file2
goto QUIT

:SKIP2
echo     Error - (%0): Second filename is missing. Usage: tfcm1 file1 file2
goto QUIT

:SKIP3
echo     Error - (%0): Third filename is present. Usage: tfcm1 file1 file2
goto QUIT

:NO1
echo     Error -(%0): File "%1" is not present.
goto QUIT

:NO2
echo     Error -(%0): File "%2" is not present.
goto QUIT

:QUIT
