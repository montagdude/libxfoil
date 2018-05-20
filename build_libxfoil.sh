#!/bin/bash

BUILDDIR=$(pwd)/build
rm -rf $BUILDDIR
mkdir $BUILDDIR

INSTALLDIR=$(pwd)/install

if [ -z "$ARCH" ]; then
  case "$( uname -m )" in
    i?86) ARCH=i586 ;;
    arm*) ARCH=arm ;;
       *) ARCH=$( uname -m ) ;;
  esac
fi

if [ "$ARCH" = "i586" ]; then
  BUILDCFLAGS="-O2 -march=i586 -mtune=i686"
  LIBDIRSUFFIX=""
elif [ "$ARCH" = "i686" ]; then
  BUILDCFLAGS="-O2 -march=i686 -mtune=i686"
  LIBDIRSUFFIX=""
elif [ "$ARCH" = "x86_64" ]; then
  BUILDCFLAGS="-O2 -fPIC"
  LIBDIRSUFFIX="64"
else
  BUILDCFLAGS="-O2"
  LIBDIRSUFFIX=""
fi

cd build
  cmake \
    -DCMAKE_C_FLAGS:STRING="$BUILDCFLAGS" \
    -DLIB_SUFFIX=${LIBDIRSUFFIX} \
    -DCMAKE_BUILD_TYPE="Release" \
    -DCMAKE_INSTALL_PREFIX=$INSTALLDIR \
    ..

  make VERBOSE=1 || exit 1
  make install || exit 1
cd ..

python setup.py build_ext
python setup.py install --prefix=$INSTALLDIR
