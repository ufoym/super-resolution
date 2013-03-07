@echo off
setLocal EnableDelayedExpansion
set name=test_hqx

if '%1'=='' goto input_error

:: compile
%~d0
cd %~dp0
call "C:\Program Files (x86)\Microsoft Visual Studio 11.0\VC\vcvarsall.bat"
cl /Fe%name% /EHsc ../image.cpp %name%.cpp
del *.obj

:: test
set syear=%date:~2,2%-%date:~5,2%-%date:~8,2%
set stime=%time:~0,2%%time:~3,2%%time:~6,2%
set dir_out="%syear% %stime% %name%"
mkdir %dir_out%
set /a i=1
for /f "tokens=*" %%f in ('dir %1 /b /a-d /s') do (
	copy %%f %dir_out%\!i!.png>nul
	%name% %%f %dir_out%\!i!.ppm>>%dir_out%\log.txt
	hqx\hqx %dir_out%\!i!.ppm %dir_out%\!i!.ppm
	set /a i+=1
	set /p="*"<nul
)

:: clean up
del %name%.exe

:: error
goto end
:input_error
echo please drag and drop your test folder...
@pause
:end