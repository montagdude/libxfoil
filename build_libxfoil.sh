#!/bin/bash

BUILDDIR=$(pwd)/build
rm -rf $BUILDDIR
mkdir $BUILDDIR

INSTALLDIR=$(pwd)/install
rm -rf $INSTALLDIR

LIBRARY_TYPE=${LIBRARY_TYPE:-"shared"}

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
  if [ "$LIBRARY_TYPE" == "shared" ]; then
    BUILDCFLAGS="-O2 -fPIC"
  else
    BUILDCFLAGS="-O2"
  fi
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
    -DLIBRARY_TYPE:STRING="$LIBRARY_TYPE" \
    ..

  make VERBOSE=1 || exit 1
  make install || exit 1
cd ..

if [ "$LIBRARY_TYPE" == "shared" ]; then
  python setup.py build_ext
  python setup.py install --prefix=$INSTALLDIR
else
  echo "Not building Python bindings - shared library is required."
fi
