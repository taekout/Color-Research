###### MAKEFILE FOR COLOR2Gray ########
########
# Variables for stuff that varies on different machines:
# Where the librarys are
# What flags to use with the C compiler
########
# comment this out for handin, this hard codes it to compile for mac

##### FIX THIS TO WORK WITH YOUR SYSTEM ############## 
INC_OPENGL =  -I/System/Library/Frameworks/Glut.framework/Versions/A/Headers/ \
		        -I/System/Library/Frameworks/OpenGL.framework/Versions/A/Headers/
LIB_OPENGL =  -framework AGL -framework OpenGL  -framework Glut
FLTK = /Users/gooch/Development/fltk-1.1.6
OTHERLIBS = -framework ApplicationServices -framework Carbon

GLDLIBS		=  -framework AGL -framework OpenGL  -framework Carbon  -framework ApplicationServices -lsupc++

POSTBUILD	= /Developer/Tools/Rez -t APPL -o
 
LDLIBS =  -L$(FLTK)/lib/ $(OTHERLIBS)   \
	  $(LIB_OPENGL) -lfltk  `fltk-config --ldflags --cxxflags --use-images` \
	  -lm -lstdc++ -lobjc 
  
CFLAGS = -g -D __APPLE__ 


########
# Some stuff that's the same for all programs and machines 
########

# The commands to call the C and C++ compilers
CC = cc 
C++ = gcc
 
# Where to find the include files:
 INCLUDE = $(INC_OPENGL) -I$(FLTK)


########
# Stuff that depends on the specific program we're compiling
# To use this file for other programs, you'll need to change things
# in this section.
########


# The name of the output program 
PROG = color2gray

# Object files that go into the final executable
OBJS = main.o image_io.o color2gray.o  gui2.o thetacircle2.o file_chooser.o

# Header files 
HDRS =  amy_colors.h images.h color2gray.h gui2.h thetacircle2.h file_chooser.h

# FLUID file rules
.fl.cxx .fl.h:	$(FLROOT)/fluid/fluid$(EXEEXT)
	echo Generating $<$(FLROOT).
	$(FLROOT)/fluid/fluid$(EXEEXT) -c $<


########
# Finally, the commands that actually make stuff
# You shouldn't have to change anything here 
# These commands are of the form:
#   output_thing:  stuff it depends on
#   <TAB> command to make output_thing 
########

 

# change a C object file whenever the corresponding source file
# or any of the header files changes
%.o : %.c $(HDRS) 
	$(CC) -c $(CFLAGS) $(INCLUDE) $<
# What this means: 
# The %.o and %.c mean any .o file and the corresponding .c file
# All .c's are recompiled if you change one of the header files
# The $< symbol refers to the FIRST dependancy - so don't try to
# put $(HDRS) before %.c!

# change a C++ object file whenever the corresponding source file
# or any of the header files changes
%.o : %.cpp $(HDRS) 
	$(C++) -c $(CFLAGS) $(INCLUDE) $<

 # change a C++ object file whenever the corresponding source file
# or any of the header files changes
%.o : %.cxx $(HDRS) 
	$(C++) -c $(CFLAGS) $(INCLUDE) $<



# re-link the program when the object files change
$(PROG):  $(OBJS)
	$(C++) $(OBJS) $(LDLIBS) -bind_at_load -o $(PROG) 


## http://www.geocities.com/rjpoling/MacOS/Makefile.txt
macosx:
	$(POSTBUILD) $(PROG) $(FLTK)/FL/mac.r