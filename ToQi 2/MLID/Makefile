# Makefile for CEM (C EM Algorithm Implementation)

CC = gcc
CFLAGS = -Wall -Wextra -O3 -std=c99 -Iinclude
LDFLAGS = -lm

# Optional GSL support (can be disabled if GSL not available)
GSL_AVAILABLE := $(shell pkg-config --exists gsl && echo yes)
ifeq ($(GSL_AVAILABLE),yes)
    CFLAGS += -DHAVE_GSL $(shell pkg-config --cflags gsl)
    LDFLAGS += $(shell pkg-config --libs gsl)
else
    $(warning GSL not found, using basic optimization only)
endif

# Optional OpenMP support
OPENMP_AVAILABLE := $(shell echo | $(CC) -fopenmp -E - 2>/dev/null && echo yes)
ifeq ($(OPENMP_AVAILABLE),yes)
    CFLAGS += -fopenmp -DHAVE_OPENMP
    LDFLAGS += -fopenmp
endif

# Debug build option
DEBUG ?= 0
ifeq ($(DEBUG),1)
    CFLAGS += -g -DDEBUG -O0
    CFLAGS := $(filter-out -O3,$(CFLAGS))
endif

# Directories
SRCDIR = src
INCDIR = include
TESTDIR = tests
BUILDDIR = build

# Source files
SOURCES = $(wildcard $(SRCDIR)/*.c)
OBJECTS = $(SOURCES:$(SRCDIR)/%.c=$(BUILDDIR)/%.o)
TARGET = mlid

# Test files
TEST_SOURCES = $(wildcard $(TESTDIR)/*.c)
TEST_OBJECTS = $(TEST_SOURCES:$(TESTDIR)/%.c=$(BUILDDIR)/test_%.o)
TEST_TARGETS = $(TEST_SOURCES:$(TESTDIR)/%.c=$(BUILDDIR)/test_%)

.PHONY: all clean test debug install

# Default target
all: $(TARGET)

# Create build directory
$(BUILDDIR):
	mkdir -p $(BUILDDIR)

# Main executable
$(TARGET): $(BUILDDIR) $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS)

# Object files
$(BUILDDIR)/%.o: $(SRCDIR)/%.c
	$(CC) $(CFLAGS) -c $< -o $@

# Test object files
$(BUILDDIR)/test_%.o: $(TESTDIR)/%.c
	$(CC) $(CFLAGS) -c $< -o $@

# Test executables
$(BUILDDIR)/test_%: $(BUILDDIR)/test_%.o $(filter-out $(BUILDDIR)/main.o, $(OBJECTS))
	$(CC) $^ -o $@ $(LDFLAGS)

# Build tests
test: $(TARGET) $(TEST_TARGETS)
	@echo "Running validation tests..."
	@for test in $(TEST_TARGETS); do \
		echo "Running $$test..."; \
		$$test || exit 1; \
	done

# Debug build
debug:
	$(MAKE) DEBUG=1

# Clean build files
clean:
	rm -rf $(BUILDDIR) $(TARGET)

# Install (optional)
install: $(TARGET)
	cp $(TARGET) /usr/local/bin/

# Show configuration
config:
	@echo "Compiler: $(CC)"
	@echo "CFLAGS: $(CFLAGS)"
	@echo "LDFLAGS: $(LDFLAGS)"
	@echo "GSL Available: $(GSL_AVAILABLE)"
	@echo "OpenMP Available: $(OPENMP_AVAILABLE)"

# Help
help:
	@echo "Available targets:"
	@echo "  all      - Build main executable (default)"
	@echo "  test     - Build and run tests"
	@echo "  debug    - Build with debug symbols"
	@echo "  clean    - Remove build files"
	@echo "  install  - Install to /usr/local/bin"
	@echo "  config   - Show build configuration"
	@echo "  help     - Show this help"
	@echo ""
	@echo "Variables:"
	@echo "  DEBUG=1  - Enable debug build"

# Dependency tracking
-include $(OBJECTS:.o=.d)
$(BUILDDIR)/%.d: $(SRCDIR)/%.c
	@$(CC) $(CFLAGS) -MM $< > $@.tmp
	@sed 's,\($*\)\.o[ :]*,$(BUILDDIR)/\1.o $@ : ,g' < $@.tmp > $@
	@rm -f $@.tmp