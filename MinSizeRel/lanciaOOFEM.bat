cd %~dp0
oofem -f %1 
REM > log.dat
REM data corrente per ultima modifica
copy /b lanciaOOFEM.bat +,,
pause