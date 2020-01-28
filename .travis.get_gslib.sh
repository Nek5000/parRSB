#!/bin/sh

VER=1.0.5

if [ -f gslib/lib/libgs.a ]; then
  echo "libgs.a found -- nothing to build."
else
  echo "Downloading gs source."
  wget https://github.com/Nek5000/gslib/archive/v${VER}.tar.gz

  echo "Building gs ..."
  tar -xvf v${VER}.tar.gz
  cd gslib-${VER}
  CC=mpicc CFLAGS="-fPIC -O2" make DESTDIR=`pwd`/../gslib install
  cd ..
  rm -rf gslib-${VER} v${VER}.tar.gz
fi
