CC=gcc
CFLAGS= -Wall -g

INCLUDE=${PWD}/include
INCLUDE_GF2=${PWD}/GF2
#SUBDIR=${PWD}/test ${PWD}/GF2 ${PWD}/GF2_extension
INCLUDE_GF2m=${PWD}/GF2_extension
INCLUDE_matrix=${PWD}/matrix
OBJECTS= test.o gf2.o gf2_test.o gf2m.o gf2m_test.o gf2_matrix.o bmatrix.o matrix_test.o
SRCS = $(OBJECTS:.o=.c)
TARGET=main.out

$(TARGET): $(OBJECTS)
	${CC} ${CFLAGS} -o $@ ${OBJECTS} -I${INCLUDE}

${OBJECTS}: test.o gf2.o gf2_test.o gf2m.o gf2m_test.o gf2_matrix.o bmatrix.o matrix_test.o
	${CC} ${CFLAGS} -c test/test.c ${INCLUDE_GF2m}/*.c ${INCLUDE_GF2}/*.c ${INCLUDE_matrix}/*.c -I${INCLUDE}

clean:
	rm -rf $(OBJECTS)