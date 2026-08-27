# Updated Makefile with proper paths
CC = gcc
CFLAGS = -O3 -fPIC -Wall -Wextra -std=c99 -march=native -flto
LDFLAGS = -shared -flto

# Source and target
SOURCE = min_cw_gbd_fq_production.c
TARGET = min_cw_gbd_fq.so

# Default target
all: $(TARGET)

# Build shared library
$(TARGET): $(SOURCE)
	@echo " Building optimized C library..."
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $<
	@echo " Built: $@ ($$(ls -lh $@ | awk '{print $$5}'))"

# Test compilation without linking
test-compile: $(SOURCE)
	@echo " Testing C compilation..."
	$(CC) $(CFLAGS) -c -o /tmp/test_compile.o $<
	@echo " Compilation successful"
	@rm -f /tmp/test_compile.o

# Debug build
debug: CFLAGS = -O0 -g -fPIC -Wall -Wextra -std=c99 -DDEBUG
debug: $(TARGET)
	@echo " Debug build complete"

# Clean build artifacts
clean:
	rm -f $(TARGET) *.o
	@echo " Cleaned build artifacts"

# Install to parent directory for easy access
install: $(TARGET)
	cp $(TARGET) ../../
	@echo " Installed $(TARGET) to project root"

# Full rebuild
rebuild: clean all

.PHONY: all test-compile debug clean install rebuild

# Build info
info:
	@echo "BUILD INFO:"
	@echo "   Source: $(SOURCE)"
	@echo "   Target: $(TARGET)"
	@echo "   Compiler: $(CC) $(CFLAGS)"
	@echo "   Features: Optimized for $(shell uname -m), LTO enabled"
