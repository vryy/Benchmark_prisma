@echo off
cls
setlocal enabledelayedexpansion
set OMP_NUM_THREADS=1
set PY_COMMAND=python

:: Ensure log directory exists
if not exist ztest_logs mkdir ztest_logs

:: Build timestamp-friendly filename (DD/MM/YYYY)
for /f "tokens=1-3 delims=/.- " %%a in ("%date%") do (
    set dd=%%a
    set mm=%%b
    set yyyy=%%c
)
set DATE_STR=%yyyy%_%mm%_%dd%

if "%~1"=="" (
    echo Run all tests
    set "OUTPUT=ztest_logs\%USERNAME%_%COMPUTERNAME%_%DATE_STR%.log"
    echo Logging to !OUTPUT!
    %PY_COMMAND% run_tests.py %PY_COMMAND% > "!OUTPUT!" 2>&1
    type "!OUTPUT!"
) else (
    echo Run tests with arguments: %*
    set "OUTPUT=ztest_logs\%USERNAME%_%COMPUTERNAME%_%DATE_STR%-%1.log"
    echo Logging to !OUTPUT!
    %PY_COMMAND% run_tests.py %PY_COMMAND% %* > "!OUTPUT!" 2>&1
    type "!OUTPUT!"
)

endlocal
