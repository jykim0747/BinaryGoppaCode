CC = gcc
LD = gcc

INCDIR = ../include
INCDIR_GF2 = ../GF2
INCDIR_GF2m = ../GF2_extension
INCDIR_matrix = ../matrix
INCDIR_McEliece = ../McEliece
INCDIR_test = ../test

VPATH = $(INCDIR) $(INCDIR_GF2) $(INCDIR_GF2m) $(INCDIR_matrix) $(INCDIR_McEliece) $(INCDIR_test)

TARGET = main.out
CFLAGS = -O2 -I$(INCDIR)

SRCS = $(foreach dir, .. $(VPATH), $(wildcard $(dir)/*.c))
SRCS := $(notdir $(SRCS))

OBJS = $(SRCS:.c=.o)

all : $(TARGET)

$(TARGET) : $(OBJS)
	$(LD) $^ -o$(TARGET)$(LIBS)
	mv $(TARGET) ../

%o:%c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -rf $(OBJS)
	rm ../$(TARGET)

depend:
	$(CC) -M $(CFLAGS) $^ > $@


