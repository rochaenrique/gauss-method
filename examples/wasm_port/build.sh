#!/bin/bash

CC="clang"
COMMONFLAGS="-Wall -Wextra -I../.."
BUILD_DIR="build"

if [[ ! -e "$BUILD_DIR" ]]
then
	mkdir $BUILD_DIR
fi

# compile basic nativly arm64
CMD="$CC $COMMONFLAGS -o $BUILD_DIR/main main.c"
echo [INFO]: Compiling 'main.c' to Darwin arm64
echo [CMD]: $CMD
eval $CMD
echo

# compile basic in wasm32
CMD="$CC $COMMONFLAGS --target=wasm32 -nostdlib -Wl,--no-entry,-export=main,-export=__heap_base -o $BUILD_DIR/main.wasm32 main.c"
echo [INFO]: Compiling 'main.c' to Wasm32
echo [CMD]: $CMD
eval $CMD
