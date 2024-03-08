@echo off
@REM make DLL(s):
nvcc -c -o cuLUT.o cuLUT.cu -DDLL
nvcc -shared -o cuLUT.dll cuLUT.o -lcudart

@REM (objdump: A good tool of gcc users to search for the addresses of functions)
@REM objdump -x cuLUT.dll

@REM make static libraries:
g++ -c  parser.cpp -o parser.o
@REM ar rcs lib_parser.a parser.o
@REM lib /OUT:lib_parser.lib parser.o 
ar rcs lib_parser.lib parser.o

@REM make main.exe
@REM the syntax of nvcc is different from gcc to link a lib
@REM g++ main.cpp -o main.exe -L. -lparser
@REM nvcc main.cu -o main.exe -L. -llib_parser -rdc=true
nvcc main.cu -o main.exe -rdc=true

@REM run
main.exe