CC = g++
CFLAGS = -std=c++11 -Wall -pedantic -g -Wno-sign-compare
EXEC_NAME = Exercice7_2022_solution
EXEC_NAME_STD = Exercice7_2022_student
INCLUDES =
LIBS =
OBJ_FILES = Exercice7_2022_solution.o
OBJ_FILES_STD = Exercice7_2022_student.o
EXEC_SUFFIX = .exe

# add suffix to exec
EXEC_NAME_SUFFIX = $(addsuffix $(EXEC_SUFFIX),$(EXEC_NAME))
EXEC_NAME_STD_SUFFIX = $(addsuffix $(EXEC_SUFFIX),$(EXEC_NAME_STD))

all : student

clean :
	rm $(EXEC_NAME_SUFFIX) $(EXEC_NAME_STD_SUFFIX)  $(OBJ_FILES) *.out
 
solution : $(OBJ_FILES)
	mkdir -p sim
	$(CC) -o $(EXEC_NAME_SUFFIX) $(OBJ_FILES) $(LIBS)

student : $(OBJ_FILES_STD)
	mkdir -p sim
	$(CC) -o $(EXEC_NAME_STD_SUFFIX) $(OBJ_FILES_STD) $(LIBS)

%.o: %.cpp
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

