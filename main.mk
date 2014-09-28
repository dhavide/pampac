#TARGET_DIR := lib
SUBMAKEFILES := src/src.mk example/KS_example.mk

# Parameters for users to configure on their systems:
CC := /usr/bin/mpicc
CFLAGS := -std=c99 -g -Wall
INCDIRS := /usr/include/gsl
