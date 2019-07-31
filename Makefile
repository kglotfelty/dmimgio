
#-----------------------------------------------------------------------
#
# 		Makefile for the tcdIO library
#
#-----------------------------------------------------------------------
MK_TOP = ../../../..
include $(MK_TOP)/Makefile.master
include $(MK_TOP)/include/Makefile.scidev

EXEC              = 
LIB_FILES         = libdmimgio.a
PAR_FILES         = 
INC_FILES         = 

SRCS	=  dmimgio.c
OBJS	= $(SRCS:.c=.o)



LIB_SRCS = $(SRCS)
LIB_OBJS = $(OBJS)

INSTALL_FILES =

include $(MK_TOP)/Makefile.all


#-----------------------------------------------------------------------
# 			MAKEFILE DEPENDENCIES	
#-----------------------------------------------------------------------

$(LIB_FILES): $(OBJS) 
	$(AR) $@ $(OBJS)
	$(RANLIB) $@
	@echo


announce1:
	@echo "   /----------------------------------------------------------\ "
	@echo "   |                 Building dmimg IO library                  | "
	@echo "   \----------------------------------------------------------/ "


