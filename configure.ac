AC_INIT(config, version-0.0.1)
AC_MSG_NOTICE([Generating Makefile for generator compilation.])

AC_LANG(C++)
AC_PROG_CXX

INC_CXX="$INC_CXX include"
SRC_CXX="$SRC_CXX src"

R_INC=`root-config --incdir`
R_LIBS=`root-config --libs --glibs` 

AC_SUBST(INCLUDE, ["$INC_CXX"])
AC_SUBST(SRC, ["$SRC_CXX"])
AC_SUBST(ROOTINC, ["$R_INC"])
AC_SUBST(ROOTLIBS, ["$R_LIBS"])

AC_OUTPUT(Makefile)
