@echo off

:: Install CHIIMP on Windows.
:: 
:: A base R install is assumed to already be present, but all dependencies
:: should be installed automatically here.

REM  Figure out path to the default R install.
for /f "delims== tokens=2 usebackq" %%x in (`ftype RWorkspace`) do set rpath=%%x
set rpath=%rpath: "%1"=%
set rpath=%rpath:"=%
set rdir=%rpath%\..\
set rscript=%rdir%\Rscript

REM  Path to chiimp source, relative to this script.
set pkgdir=%~dp0

REM  Run bulk of the install within R.
"%rscript%" --vanilla "%pkgdir%\inst\installer\install.R"
pause
