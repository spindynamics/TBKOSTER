# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

# Default target executed when no arguments are given to make.
default_target: all

.PHONY : default_target

# Allow only one "make -f Makefile2" at a time, but pass parallelism.
.NOTPARALLEL:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /home/barreto/bin/cmake-3.10.2-Linux-x86_64/bin/cmake

# The command to remove a file.
RM = /home/barreto/bin/cmake-3.10.2-Linux-x86_64/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/barreto/github/DyNaMol

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/barreto/github/DyNaMol

#=============================================================================
# Targets provided globally by CMake.

# Special rule for the target edit_cache
edit_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake cache editor..."
	/home/barreto/bin/cmake-3.10.2-Linux-x86_64/bin/ccmake -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : edit_cache

# Special rule for the target edit_cache
edit_cache/fast: edit_cache

.PHONY : edit_cache/fast

# Special rule for the target rebuild_cache
rebuild_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake to regenerate build system..."
	/home/barreto/bin/cmake-3.10.2-Linux-x86_64/bin/cmake -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : rebuild_cache

# Special rule for the target rebuild_cache
rebuild_cache/fast: rebuild_cache

.PHONY : rebuild_cache/fast

# The main all target
all: cmake_check_build_system
	$(CMAKE_COMMAND) -E cmake_progress_start /home/barreto/github/DyNaMol/CMakeFiles /home/barreto/github/DyNaMol/CMakeFiles/progress.marks
	$(MAKE) -f CMakeFiles/Makefile2 all
	$(CMAKE_COMMAND) -E cmake_progress_start /home/barreto/github/DyNaMol/CMakeFiles 0
.PHONY : all

# The main clean target
clean:
	$(MAKE) -f CMakeFiles/Makefile2 clean
.PHONY : clean

# The main clean target
clean/fast: clean

.PHONY : clean/fast

# Prepare targets for installation.
preinstall: all
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall

# Prepare targets for installation.
preinstall/fast:
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall/fast

# clear depends
depend:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 1
.PHONY : depend

#=============================================================================
# Target rules for targets named doc

# Build rule for target.
doc: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 doc
.PHONY : doc

# fast build rule for target.
doc/fast:
	$(MAKE) -f CMakeFiles/doc.dir/build.make CMakeFiles/doc.dir/build
.PHONY : doc/fast

#=============================================================================
# Target rules for targets named UserManual_auxclean

# Build rule for target.
UserManual_auxclean: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 UserManual_auxclean
.PHONY : UserManual_auxclean

# fast build rule for target.
UserManual_auxclean/fast:
	$(MAKE) -f CMakeFiles/UserManual_auxclean.dir/build.make CMakeFiles/UserManual_auxclean.dir/build
.PHONY : UserManual_auxclean/fast

#=============================================================================
# Target rules for targets named _UserManual

# Build rule for target.
_UserManual: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 _UserManual
.PHONY : _UserManual

# fast build rule for target.
_UserManual/fast:
	$(MAKE) -f CMakeFiles/_UserManual.dir/build.make CMakeFiles/_UserManual.dir/build
.PHONY : _UserManual/fast

#=============================================================================
# Target rules for targets named UserManual_safepdf

# Build rule for target.
UserManual_safepdf: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 UserManual_safepdf
.PHONY : UserManual_safepdf

# fast build rule for target.
UserManual_safepdf/fast:
	$(MAKE) -f CMakeFiles/UserManual_safepdf.dir/build.make CMakeFiles/UserManual_safepdf.dir/build
.PHONY : UserManual_safepdf/fast

#=============================================================================
# Target rules for targets named UserManual

# Build rule for target.
UserManual: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 UserManual
.PHONY : UserManual

# fast build rule for target.
UserManual/fast:
	$(MAKE) -f CMakeFiles/UserManual.dir/build.make CMakeFiles/UserManual.dir/build
.PHONY : UserManual/fast

#=============================================================================
# Target rules for targets named dvi

# Build rule for target.
dvi: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 dvi
.PHONY : dvi

# fast build rule for target.
dvi/fast:
	$(MAKE) -f CMakeFiles/dvi.dir/build.make CMakeFiles/dvi.dir/build
.PHONY : dvi/fast

#=============================================================================
# Target rules for targets named ps

# Build rule for target.
ps: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 ps
.PHONY : ps

# fast build rule for target.
ps/fast:
	$(MAKE) -f CMakeFiles/ps.dir/build.make CMakeFiles/ps.dir/build
.PHONY : ps/fast

#=============================================================================
# Target rules for targets named distclean

# Build rule for target.
distclean: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 distclean
.PHONY : distclean

# fast build rule for target.
distclean/fast:
	$(MAKE) -f CMakeFiles/distclean.dir/build.make CMakeFiles/distclean.dir/build
.PHONY : distclean/fast

#=============================================================================
# Target rules for targets named safepdf

# Build rule for target.
safepdf: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 safepdf
.PHONY : safepdf

# fast build rule for target.
safepdf/fast:
	$(MAKE) -f CMakeFiles/safepdf.dir/build.make CMakeFiles/safepdf.dir/build
.PHONY : safepdf/fast

#=============================================================================
# Target rules for targets named pdf

# Build rule for target.
pdf: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 pdf
.PHONY : pdf

# fast build rule for target.
pdf/fast:
	$(MAKE) -f CMakeFiles/pdf.dir/build.make CMakeFiles/pdf.dir/build
.PHONY : pdf/fast

#=============================================================================
# Target rules for targets named UserManual_pdf

# Build rule for target.
UserManual_pdf: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 UserManual_pdf
.PHONY : UserManual_pdf

# fast build rule for target.
UserManual_pdf/fast:
	$(MAKE) -f CMakeFiles/UserManual_pdf.dir/build.make CMakeFiles/UserManual_pdf.dir/build
.PHONY : UserManual_pdf/fast

#=============================================================================
# Target rules for targets named UserManual_dvi

# Build rule for target.
UserManual_dvi: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 UserManual_dvi
.PHONY : UserManual_dvi

# fast build rule for target.
UserManual_dvi/fast:
	$(MAKE) -f CMakeFiles/UserManual_dvi.dir/build.make CMakeFiles/UserManual_dvi.dir/build
.PHONY : UserManual_dvi/fast

#=============================================================================
# Target rules for targets named html

# Build rule for target.
html: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 html
.PHONY : html

# fast build rule for target.
html/fast:
	$(MAKE) -f CMakeFiles/html.dir/build.make CMakeFiles/html.dir/build
.PHONY : html/fast

#=============================================================================
# Target rules for targets named auxclean

# Build rule for target.
auxclean: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 auxclean
.PHONY : auxclean

# fast build rule for target.
auxclean/fast:
	$(MAKE) -f CMakeFiles/auxclean.dir/build.make CMakeFiles/auxclean.dir/build
.PHONY : auxclean/fast

#=============================================================================
# Target rules for targets named UserManual_ps

# Build rule for target.
UserManual_ps: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 UserManual_ps
.PHONY : UserManual_ps

# fast build rule for target.
UserManual_ps/fast:
	$(MAKE) -f CMakeFiles/UserManual_ps.dir/build.make CMakeFiles/UserManual_ps.dir/build
.PHONY : UserManual_ps/fast

#=============================================================================
# Target rules for targets named DyNaMol

# Build rule for target.
DyNaMol: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 DyNaMol
.PHONY : DyNaMol

# fast build rule for target.
DyNaMol/fast:
	$(MAKE) -f src/CMakeFiles/DyNaMol.dir/build.make src/CMakeFiles/DyNaMol.dir/build
.PHONY : DyNaMol/fast

#=============================================================================
# Target rules for targets named pdos.x

# Build rule for target.
pdos.x: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 pdos.x
.PHONY : pdos.x

# fast build rule for target.
pdos.x/fast:
	$(MAKE) -f src/CMakeFiles/pdos.x.dir/build.make src/CMakeFiles/pdos.x.dir/build
.PHONY : pdos.x/fast

#=============================================================================
# Target rules for targets named DyNaMol.x

# Build rule for target.
DyNaMol.x: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 DyNaMol.x
.PHONY : DyNaMol.x

# fast build rule for target.
DyNaMol.x/fast:
	$(MAKE) -f src/CMakeFiles/DyNaMol.x.dir/build.make src/CMakeFiles/DyNaMol.x.dir/build
.PHONY : DyNaMol.x/fast

#=============================================================================
# Target rules for targets named bands.x

# Build rule for target.
bands.x: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 bands.x
.PHONY : bands.x

# fast build rule for target.
bands.x/fast:
	$(MAKE) -f src/CMakeFiles/bands.x.dir/build.make src/CMakeFiles/bands.x.dir/build
.PHONY : bands.x/fast

# Help Target
help:
	@echo "The following are some of the valid targets for this Makefile:"
	@echo "... all (the default if no target is provided)"
	@echo "... clean"
	@echo "... depend"
	@echo "... edit_cache"
	@echo "... doc"
	@echo "... UserManual_auxclean"
	@echo "... _UserManual"
	@echo "... UserManual_safepdf"
	@echo "... UserManual"
	@echo "... dvi"
	@echo "... ps"
	@echo "... rebuild_cache"
	@echo "... distclean"
	@echo "... safepdf"
	@echo "... pdf"
	@echo "... UserManual_pdf"
	@echo "... UserManual_dvi"
	@echo "... html"
	@echo "... auxclean"
	@echo "... UserManual_ps"
	@echo "... DyNaMol"
	@echo "... pdos.x"
	@echo "... DyNaMol.x"
	@echo "... bands.x"
.PHONY : help



#=============================================================================
# Special targets to cleanup operation of make.

# Special rule to run CMake to check the build system integrity.
# No rule that depends on this can have commands that come from listfiles
# because they might be regenerated.
cmake_check_build_system:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 0
.PHONY : cmake_check_build_system

