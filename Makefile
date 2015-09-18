 CPPC = g++
 FLAGS = -O3

OBJ = fcyl2crt.o

fcyl2crt: $(OBJ)
	$(CPPC) $(LFLAGS) -o fcyl2crt $(OBJ)

clean:
	rm -f *.o fcyl2crt

%.o: %.cpp
	$(CPPC) $(FLAGS) -c $^ -o $@
