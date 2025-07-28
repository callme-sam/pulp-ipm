CC = gcc
CFLAGS = -Wall -Wextra -Iinc
LDFLAGS = -lgsl -lgslcblas -lm

SRC = main.c $(wildcard src/*.c)
OBJ = $(patsubst %.c, build/obj/%.o, $(SRC))
TARGET = build/ipm_solver

# Variabile per verbose (default off)
VERBOSE ?= 0

ifeq ($(VERBOSE),1)
  V =
else
  V = @
endif

# Crea le directory necessarie
$(shell mkdir -p build/obj/src)

all: $(TARGET)

# Regola per l'eseguibile
$(TARGET): $(OBJ)
	$(V)$(CC) -o $@ $^ $(LDFLAGS)
	$(if $(filter 1,$(VERBOSE)), echo "Linking $@")

# Compila ogni .c in build/obj
build/obj/%.o: %.c inc/ipm.h
	@mkdir -p $(dir $@)
	$(V)$(CC) $(CFLAGS) -c $< -o $@
	$(if $(filter 1,$(VERBOSE)), echo "Compiling $<")

clean:
	$(V)rm -rf build/
	$(if $(filter 1,$(VERBOSE)), echo "Cleaned build directory")

run: $(TARGET)
	$(V)./$(TARGET)
