@echo off
rem tchs.bat

rem This is a DOS batch file to run an H2OI95 test
rem case. It is used by tch.bat. This batch file in turn
rem uses tfcm1.bat to compare output files. This particular
rem comparison disregards the first line, which contains
rem the H2OI95 build information.

rem Last revised 05/18/2020 by TJW

rem Note: in setting codeDIR, it may be necessary to change
rem the drive letter.

rem set codeDIR=C:\H2OI95_v1.1
rem set codeEXE=%codeDIR%\Bin\H2OI95.exe

rem The variables codeDIR and codeEXE should have been set
rem in tch.bat.

if exist %codeEXE% goto OK1
  echo     Error -(%0): Code executable "%codeEXE%" is not present.
  goto QUIT
:OK1

rem The variable batfs2 should have been set in tch.bat.

if exist %batfs2% goto OK2
  echo   Error - (%0): Batch file %batfs2% is missing.
  goto QUIT
:OK2

rem Here %1 is the name of the specified test case
rem folder, which corresponds a single test case.

if %1.==. goto SKIP1

rem if NOT exist %1 goto SKIP2

echo Processing test case %1

rem echo The current folder is:
rem pwd

rem Change to the test case folder.
rem echo Change to the test case folder.
cd %1
rem echo The current folder is:
rem pwd

rem A valid test case folder must contain a file
rem named input, or an archival copy thereof.

rem Get the test case string without the leading
rem two characters, which denote the test case type.

set str=%1
set str=%str:~2,5%
echo   The short string is %str%


set inputa=input_%str%
rem echo   The archival input filename is %inputa%

if NOT exist input goto NOINP
  echo   The normal input file is present.

rem XXXXX
rem touch input
rem XXXXX

  set exinp=TRUE
  goto END00
:NOINP
  echo   The normal input file is not present.
  set exinp=FALSE
:END00

if NOT exist %inputa% goto NOINPA
  echo   The archival input file %inputa% is present.

rem XXXXX
rem touch %inputa%
rem XXXXX

  set exinpa=TRUE
  goto END01
:NOINPA
  echo   The archival input file %inputa% is not present.
  set exinpa=FALSE
:END01

if %exinp%==FALSE goto END05

rem The normal input file is present.

  if %exinpa%==FALSE goto END04

rem The archival input file is also present.
rem Compare the normal and archival input files.

    fc input %inputa% > nul
    if %errorlevel% == 0 goto OKAY0
      echo   The files input and %inputa% are different
      goto ENDOKAY0
    :OKAY0
      echo   The files input and %inputa% are identical
    :ENDOKAY0
    goto INP2F
  :END04
:END05

rem The normal input file is not present.

  if %exinpa%==FALSE goto NOINRES

rem The archival input file is present. Copy it to
rem the normal input file.

    echo   Copying %inputa% to input
    copy %inputa% input
    set exinp=TRUE
    goto INP2F
  :NOINRES
:INP2F

IF %exinp%==TRUE goto RUN0
   echo   An input file could not be found
   goto NOTEST
:RUN0


set outputa=output_%str%
set outa=out_%str%
set mtaba=mtab_%str%.csv
set ctaba=ctab_%str%.csv
set xtaba=xtab_%str%.csv

rem The archival output file is %outputa%

set exoutpa=FALSE
if NOT exist %outputa% goto NOOUTPA
  echo   The archival output file %outputa% is present.
  set exoutpa=TRUE
  goto END07
:NOOUTPA
  echo   The archival output file %outputa% is not present.
  echo   No test comparison will be made.
:END07

rem The archival out filename is %outa%

set exouta=FALSE
if NOT exist %outa% goto NOOUTA
  echo   The archival out file %outa% is present.
  set exouta=TRUE
  goto END09
:NOOUTA
  echo   The archival out file %outa% is not present.
  echo   No test comparison will be made.
:END09

rem The archival mtab filename is %mtaba%

set exmtba=FALSE
if NOT exist %mtaba% goto NOMTBA
  echo   The archival mtab file %mtaba% is present.
  set exmtba=TRUE
  goto END10
:NOMTBA
  echo   The archival mtab file %mtaba% is not present.
  echo   No test comparison will be made.
:END10

rem The archival ctab filename is %ctaba%

set exctba=FALSE
if NOT exist %ctaba% goto NOCTBA
  echo   The archival ctab file %ctaba% is present.
  set exctba=TRUE
  goto END11
:NOCTBA
  echo   The archival ctab file %ctaba% is not present.
  echo   No test comparison will be made.
:END11

rem The archival xtab filename is %xtaba%

set exxtba=FALSE
if NOT exist %xtaba% goto NOXTBA
  echo   The archival xtab file %xtaba% is present.
  set exxtba=TRUE
  goto END12
:NOXTBA
  echo   The archival xtab file %xtaba% is not present.
  echo   No test comparison will be made.
:END12

echo   -----

rem Delete any existing normal output files.

if NOT exist output goto NOOUTP
  echo   Deleting existing output file
  del output
:NOOUTP

if NOT exist out goto NOOUT
  echo   Deleting existing out file
  del out
:NOOUT

if NOT exist mtab.csv goto END20
  echo   Deleting existing mtab.csv file.
  del mtab.csv
:END20

if NOT exist ctab.csv goto END30
  echo   Deleting existing ctab.csv file.
  del ctab.csv
:END30

if NOT exist xtab.csv goto END40
  echo   Deleting existing xtab.csv file.
  del xtab.csv
:END40

echo   -----

echo   Running H2OI95
%codeEXE% > out

echo   -----

rem Compare normal and archival "output" files.

set errstr_grp=FALSE

rem Compare the normal and archival output files

if NOT exist %outputa% goto NO1
  CALL %batfs2% output %outputa%
  if %errstr% == FALSE goto OKAY1
    set errstr_grp=TRUE
    echo   --------------------
    fc output %outputa%
rem diff output %outputa%
    echo   --------------------
  :OKAY1
:NO1

rem Compare the normal and archival out files

if NOT exist %outa% goto NO2
  CALL %batfs2% out %outa%
  if %errstr% == FALSE goto OKAY2
    set errstr_grp=TRUE
    echo   --------------------
    fc out %outa%
rem diff out %outa%
    echo   --------------------
  :OKAY2
:NO2

rem Compare the normal and archival mtab files

if NOT exist %mtaba% goto NO3
  CALL %batfs2% mtab.csv %mtaba%
  if %errstr% == FALSE goto OKAY3
    set errstr_grp=TRUE
    echo   --------------------
    fc mtab.csv %mtaba%
rem diff mtab.csv %mtaba%
    echo   --------------------
  :OKAY3
:NO3

rem Compare the normal and archival ctab files

if NOT exist %ctaba% goto NO4
  CALL %batfs2% ctab.csv %ctaba%
  if %errstr% == FALSE goto OKAY4
    set errstr_grp=TRUE
    echo   --------------------
    fc ctab.csv %ctaba%
rem diff ctab.csv %ctaba%
    echo   --------------------
  :OKAY4
:NO4

rem Compare the normal and archival xtab files

if NOT exist %xtaba% goto NO5
  CALL %batfs2% xtab.csv %xtaba%
  if %errstr% == FALSE goto OKAY5
    set errstr_grp=TRUE
    echo   --------------------
    fc xtab.csv %xtaba%
rem diff xtab.csv %xtaba%
    echo   --------------------
  :OKAY5
:NO5

if %errstr_grp%==FALSE goto ENDOK8
echo   ----- DIFFERENCES DETECTED -----
goto ENDOK9
:ENDOK8
echo   NO DIFFERENCES FOUND
:ENDOK9



rem The following line if active skips updating archival output files.
rem This is the normal testing procedure (not updating).

rem XXXXX
goto SKIP6
rem XXXXX

echo Update archival output files

echo del %outputa%
del %outputa%
echo ren output %outputa%
ren output %outputa%

echo del %outa%
del %outa%
echo ren out %outa%
ren out %outa%

echo del %mtaba%
del %mtaba%
echo ren mtab.csv %mtaba%
ren mtab.csv %mtaba%

echo del %ctaba%
del %ctaba%
echo ren ctab.csv %ctaba%
ren ctab.csv %ctaba%

echo del %xtaba%
del %xtaba%
echo ren xtab.csv %xtaba%
ren xtab.csv %xtaba%

rem If archival files have been updated, the normal output files
rem no longer exist as such.

goto SKIP7

:SKIP6



rem The following line if active skips removing the normal output files.
rem This is the normal testing procedure (not removing).

rem XXXXX
goto SKIP7
rem XXXXX

echo Remove the normal output files.

echo del output
del output
echo del out
del out

echo del mtab.csv
del mtab.csv
echo del ctab.csv
del ctab.csv
echo del xtab.csv
del xtab.csv

:SKIP7



goto QUIT

:SKIP1
echo     Error - (%0): Missing argument.
goto QUIT

:SKIP2
echo     Error -(%0): Folder "%1" is not present.
goto QUIT

:NOTEST
echo   This may not be a test case folder.
goto QUIT

:QUIT

rem Return to the parent folder.
cd ..

rem End of this batch file
