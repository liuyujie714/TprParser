@echo off
call "C:\Users\liuyujie714\Desktop\emsdk\emsdk_env.bat"

rem emcc -lembind main.cpp Reader.cpp Utils.cpp -o test.html -sNO_DISABLE_EXCEPTION_CATCHING --preload-file .\test\benchBFC_FEP.tpr && emrun test.html

emcc -lembind -O2 Reader.cpp Utils.cpp emcc.cpp -o libgmxtpr.js -sNO_DISABLE_EXCEPTION_CATCHING -sFORCE_FILESYSTEM -sMODULARIZE=1 -sEXPORT_NAME="TprModule" -sEXPORTED_RUNTIME_METHODS=['FS'] -lworkerfs.js -sWASM=0 -sALLOW_MEMORY_GROWTH=1

