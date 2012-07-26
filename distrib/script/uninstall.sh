#!/bin/sh

echo uninstalling ibex...

IBEXPATH=`pwd`

TARGET_DIR=/usr/local/include
TARGET_LINK=$TARGET_DIR/ibex
SOURCE_LINK=$IBEXPATH/include

# ============================================== 
./script/unlink.sh $SOURCE_LINK $TARGET_DIR $TARGET_LINK 

if [ $? != 0 ]; then
    exit
fi
# ============================================== 

TARGET_DIR=/usr/local/lib
TARGET_LINK=$TARGET_DIR/libibex.so
SOURCE_LINK=$IBEXPATH/lib/libibex.so

# ============================================== 
./script/unlink.sh $SOURCE_LINK $TARGET_DIR $TARGET_LINK 

if [ $? != 0 ]; then
    exit
fi
# ============================================== 

echo -n \\t Running ldconfig ..........
ldconfig
echo done.

echo done.