#!/bin/bash

if [[ ! -e "build" ]]
then
	mkdir build
fi
   
COMMONFLAGS="-Wall -Wextra -I."
CC="clang"

# compile basic nativly arm64
CMD="$CC $COMMONFLAGS -o build/basic examples/basic.c"
echo [INFO]: Compiling 'example/basic.c' to Darwin arm64
echo [CMD]: $eval
CMD $CMD

# compile basic in wasm
