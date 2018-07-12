AC_INIT([generator], [version-1.0], [jhoskins@jlab.org])
AC_MSG_NOTICE([Generating Makefile for generator compilation.])

AC_LANG(C++)
AC_PROG_CXX

# introduce the optional configure parameter for a non-standard install prefix of XXX
AC_ARG_WITH([root],
    [AS_HELP_STRING([--with-root],
        [Compile with root support.])],
    [with_root=1],
    [with_root=0])
    
# if optional parameter used, extend path flags for compliler and linker
if test $with_root = 1 ; then
    # extend the compiler and linker flags according to the path set
    AC_CXXFLAGS="$AC_CXXFLAGS `root-config --incdir`"
    AC_LDFLAGS="$AC_LDFLAGS `root-config --libs --glibs`"
fi

AC_SUBST([AC_CXXFLAGS])
AC_SUBST([AC_LDFLAGS])

AC_OUTPUT(Makefile)
