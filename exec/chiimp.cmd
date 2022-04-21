@echo off

:: Windows wrapper to the CHIIMP command-line script.
::
:: This will locate the R and RStudio installations to use for CHIIMP.  This
:: script will wait for a keypress before exiting since it presumably opened
:: its own terminal window.

REM  This directory
set dir=%~dp0

REM  Figure out path to the default R install.
for /f "delims== tokens=2 usebackq" %%x in (`ftype RWorkspace`) do set rpath=%%x
set rpath=%rpath: "%1"=%
set rpath=%rpath:"=%
set rdir=%rpath%\..\

for /f "tokens=* usebackq" %%x in (`%rdir%\Rscript %dir%\find_pandoc.R`) do set RSTUDIO_PANDOC=%%x

if "%~1"=="" (
	echo.To run CHIIMP, drag and drop a configuration file onto this icon.
	echo.
	echo.For more information see the user guide bundled with the program or here:
	echo.https://shawhahnlab.github.io/chiimp/GUIDE.pdf
	echo.
) else (
	"%rdir%\RScript" "%dir%\chiimp" %*
)
pause
