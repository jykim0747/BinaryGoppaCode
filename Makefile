CC=gcc
CFLAGS= -Wall -g

INCLUDE=${PWD}/include
INCLUDE_GF2=${PWD}/GF2
INCLUDE_GF2m=${PWD}/GF2_extension
OBJECTS= test.o gf2.o gf2_test.o gf2m.o gf2m_test.o
SRCS = $(OBJECTS:.o=.c)
TARGET=main.out


$(TARGET): $(OBJECTS)
	${CC} ${CFLAGS} -o $@ ${OBJECTS} -I${INCLUDE}

${OBJECTS}: test.o gf2.o gf2_test.o gf2m.o gf2m_test.o
	${CC} ${CFLAGS} -c test/test.c ${INCLUDE_GF2m}/*.c ${INCLUDE_GF2}/*.c -I${INCLUDE}

clean:
	rm -rf $(OBJECTS)