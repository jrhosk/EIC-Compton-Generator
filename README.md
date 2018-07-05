# EIC-Compton-Generator

# Installation

The generator sotware now uses autoconf to build and compile the code. 

1) Initialize autoconf:
>>>> autoconf
2) Run the configuration 
>>>> ./configure
     
     The configuration will look at the template, Makefile.in, to build a makefile. Changes to the form of the makefile 
     should be done here. Changes to the configuration should be done in configure.ac. I don't advise making changes 
     unless you are very confident you know what you are doing.

     You can overwrite options with in the Makefile using command-line option with configuration. For instance if you
     wanted to compile the code using c++11 (no currently supported) you would use the CPPFLAGS.

>>>> ./configure CXX='g++ -std=c++11'
     
     Compilers flags can be change in a similar fashion using CXXFLAGS="".

2a) Compiling with options (ROOT)
    
    Currently the only option here is to compile including ROOT support, which is required until the next release. 
    Compilation including root is achieved using the --with-root flag.

>>>> ./configure --with-root

3) If configure finished with no errors. The code can be compiled using

>>>> make

     The generator should noe me ready to use. The generator options can be seen using.

>>>> ./generator --help

     Plesae direct any problems to myself in the form of a bug report or via email at jhoskins@jlab.org.
