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
set rexe=%rdir%\R

REM  Path to chiimp source, relative to this script.
set pkgdir=%~dp0
REM  Since I'm pasting text directly into R commands below we'll need to 
REM  escape the backslashes.
REM  Also, adding a stub to the end of the string to work around Windows' 
REM  weirdness to do with trailing slashes as described here:
REM  https://github.com/hadley/devtools/issues/1614
set pkgdir_r=%pkgdir:\=\\%\\.

REM  R commands to run below.  Trying to avoid confusing cmd.exe by leaving
REM  out spaces and using single-quotes throughout.  If this script grows any
REM  more complex it might be worth switching to PowerShell.
set devtools_setup=install.packages('devtools',repos='https://cloud.r-project.org')
set bioclite_setup=source('https://bioconductor.org/biocLite.R');biocLite();biocLite('msa')
set deps_setup=devtools::install_deps('%pkgdir_r%',dependencies=TRUE)
set chiimp_test=quit(save='no',status=sum(as.data.frame(devtools::test('%pkgdir_r%'))$failed))
set chiimp_setup=devtools::install('%pkgdir_r%')
set chiimp_get_path=cat(system.file('bin','chiimp',package='chiimp'))

"%rexe%" --version
echo.
echo ### Installing devtools
echo.
"%rexe%" --slave -e %devtools_setup%
echo.
echo ### Installing Bioconductor and MSA
echo.
"%rexe%" --slave -e %bioclite_setup%
echo.
echo ### Installing dependencies
echo.
"%rexe%" --slave -e %deps_setup%
echo.
echo ### Testing CHIIMP
echo.
"%rexe%" --slave -e %chiimp_test%
if errorlevel 1 (
	echo.
	echo.
	echo     Warning: Tests indicated failures.
	echo.
	echo.
)
echo.
echo ### Installing CHIIMP
echo.
"%rexe%" --slave -e %chiimp_setup%

for /f %%x in ('"%rexe%" --slave -e cat^(system.file^('bin'^,'chiimp.cmd'^,package^='chiimp'^)^)') do set chiimp_path=%%x
REM  https://stackoverflow.com/a/30029955/6073858
powershell "$s=(New-Object -COM WScript.Shell).CreateShortcut('%userprofile%\Desktop\CHIIMP.lnk');$s.TargetPath='%chiimp_path%';$s.Save()"
