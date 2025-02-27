
CC = gcc
CFLAGS = -Wall -O2 -std=c99
LDFLAGS = -lm -lopenblas


SRC = main.c matrix_decomposition_dynamic.c linear_sys_equs.c matrix_funcs.c ODE_solver.c operators.c
OBJ = $(SRC:.c=.o)
TARGET = my_program


all: $(TARGET)


$(TARGET): $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)


%.o: %.c
	$(CC) $(CFLAGS) -c $<


clean:
	rm -f $(OBJ) $(TARGET)

.PHONY: all clean
