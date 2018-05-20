SET INSTALLDIR=%CD%\windows
IF EXIST %INSTALLDIR% RMDIR /S /Q %INSTALLDIR%
IF EXIST build RMDIR /S /Q build
MKDIR %INSTALLDIR%
MKDIR build
SET BUILDCFLAGS="-O2 -fPIC"

CD build
cmake -G "MinGW Makefiles" ^
  -DCMAKE_C_FLAGS:STRING=%BUILDCFLAGS% ^
  -DCMAKE_BUILD_TYPE:STRING="Release" ^
  -DCMAKE_INSTALL_PREFIX=%INSTALLDIR% ^
  ..

mingw32-make VERBOSE=1
mingw32-make install
CD ..

python setup.py build_ext
python setup.py install --prefix=$INSTALLDIR
