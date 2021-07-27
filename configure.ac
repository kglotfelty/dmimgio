AC_INIT([dmimgio], [4.13.5], [glotfeltyk@si.edu])
AM_INIT_AUTOMAKE([foreign -Wall])
AM_MAINTAINER_MODE([disable])
AC_PROG_CC
AC_PROG_RANLIB
AM_PROG_CC_C_O
AM_PROG_AR

# -- Hack to provide missing pkg-config files


export PKG_CONFIG_PATH=$prefix/lib/pkgconfig:$prefix/ots/lib/pkgconfig:$PKG_CONFIG_PATH
PKG_CONFIG_PATH=./lib/pkgconfig:${PKG_CONFIG_PATH}

PKG_CHECK_MODULES( [CIAO], [
  ascdm
  cxcparam
  region
  err
  ds
])

#AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([
  Makefile
  src/Makefile
  src/dmimgio.pc
])


AC_OUTPUT
  
