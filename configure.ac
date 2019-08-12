AC_INIT([dmimgio], [1.0.0], [glotfeltyk@si.edu])
AM_INIT_AUTOMAKE([foreign -Wall])
AC_PROG_CC
AC_PROG_RANLIB
AM_PROG_CC_C_O
AM_PROG_AR

# -- Hack to provide missing pkg-config files

sed s,@prefix@,$prefix,g < lib/pkgconfig/ascdm.pc.in > lib/pkgconfig/ascdm.pc
sed s,@prefix@,$prefix,g < lib/pkgconfig/readline.pc.in > lib/pkgconfig/readline.pc
sed s,@prefix@,$prefix,g < lib/pkgconfig/region.pc.in > lib/pkgconfig/region.pc


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
  
