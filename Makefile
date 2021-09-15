CC = gcc
CFLAGS = -Wall -g

OBJDIR = ./obj

all:
	cd $(OBJDIR) && make

clean:
	cd $(OBJDIR) && make clean

depend:
	cd $(OBJDIR) && make depend