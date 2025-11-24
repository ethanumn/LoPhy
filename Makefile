MAKE_PATH=$(realpath $(shell dirname $(firstword $(MAKEFILE_LIST))))
CPPFLAGS=-I$(MAKE_PATH)/include 
SRC_DIR=$(MAKE_PATH)/src

CC=g++ -std=c++20 -O3 -Wall -L/opt/homebrew/lib -Wextra -pedantic -lm -Wl,-rpath,/opt/homebrew/lib
$(info $(MAKE_PATH))

LoPhy: $(wildcard $(SRC_DIR)/*.cpp)
	    $(CC) $(CPPFLAGS) -o $(MAKE_PATH)/bin/LoPhy $(wildcard $(SRC_DIR)/*.cpp)
