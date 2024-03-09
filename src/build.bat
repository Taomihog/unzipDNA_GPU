@echo off
@echo Cleaning up directory:
set "batch_folder=%~dp0"
cd %batch_folder%
del *.exp
del *.lib
del *.o
del *.dll
del *.exe
del *.csv

@echo Making cuLUT.dll:
@REM I am not sure about what the "cudart" do
nvcc -c -o cuLUT.o cuLUT.cu -DDLL
nvcc -shared -o cuLUT.dll cuLUT.o -lcudart

@echo Making cuUnzip.dll:
nvcc -c -o cuUnzip.o cuUnzip.cu -DDLL
nvcc -shared -o cuUnzip.dll cuUnzip.o -lcudart

@echo (Note: objdump -- A good tool of gcc users to search for the addresses of functions)
@REM objdump -x cuLUT.dll

@echo Making main:
@echo (Note: the syntax of nvcc is different from gcc to link a lib)
@REM g++ main.cpp -o main.exe -L. -lparser
@REM nvcc main.cu -o main.exe -L. -llib_parser -rdc=true
nvcc main.cu -o main.exe -rdc=true -std=c++17

@echo Cleaning up directory:
del *.exp
del *.lib
del *.o

@echo Testing: 
IF EXIST "../examples/parsed" (
    rd /s /q "../examples/parsed_unzip_curves"
    main ../examples/parsed 
)