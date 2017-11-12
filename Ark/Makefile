vpath %.h Header:Lib
vpath %.cpp Cpp:Lib

QModel = D2V16

CC = g++

CFLAGS = -std=c++11 -I./Header -I./Lib -O3

MPFLAG = -fopenmp

DIR_SRC = ./Cpp

DIR_OBJ = ./Obj

SRC = $(wildcard ${DIR_SRC}/*.cpp)

objects = $(patsubst %.cpp,${DIR_OBJ}/%.o,$(notdir ${SRC})) ${DIR_OBJ}/${QModel}.o

main:$(objects)

	$(CC) $(objects) -o $@ ${MPFLAG}

${DIR_OBJ}/%.o:${DIR_SRC}/%.cpp

	$(CC) $(CFLAGS) -c $< -o $@ ${MPFLAG}

${DIR_OBJ}/${QModel}.o:./Lib/${QModel}.cpp

	$(CC) $(CFLAGS) -c $< -o $@

.PHONY:clean
clean:
	find ${DIR_OBJ} -name "*.o" -exec mv {} ~/.trash/ \;

# headers = ZeroDeclaration.h Mesh_2D.h ZeroConstant.h $(QModel).h ZeroReference.h ZeroFlip.h

# main : $(objects)
# 	$(CC) -o main $(objects) -fopenmp

# main.o : main.cpp $(headers)
# 	$(CC) -c -O3 -std=c++11 main.cpp

# MeshConstruct.o : MeshConstruct.cpp MeshConstructFunc.h $(headers)
# 	$(CC) -c -O3 -std=c++11 MeshConstruct.cpp

# Preprocess.o : Preprocess.cpp $(headers)
# 	$(CC) -c -O3 -std=c++11 Preprocess.cpp

# DUGKS2DSolver.o : DUGKS2DSolver.cpp NewtonCotes.h $(headers)
# 	$(CC) -c -O3 -std=c++11 DUGKS2DSolver.cpp -fopenmp

# $(QModel).o : $(QModel).cpp $(headers)
# 	$(CC) -c -O3 -std=c++11 $(QModel).cpp

# Mesh_2D.o : Mesh_2D.cpp $(headers)
# 	$(CC) -c -O3 -std=c++11 Mesh_2D.cpp

# Output.o : Output.cpp $(headers)
# 	$(CC) -c -O3 -std=c++11 Output.cpp

# .PHONY : clean
# clean :
# 	rm pexe $(objects)