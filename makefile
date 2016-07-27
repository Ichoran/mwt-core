# Copyright Nicholas A. Swierczek, Rex A. Kerr and HHMI Janelia, 2007-2015.
# Copyright Rex A. Kerr and Calico Life Sciences, 2015-2016.
# This file is covered by the LGPL 2.1 license.

SHELL=/bin/bash

ifeq ($(OS),Windows_NT)
    CC = mingw32-g++ -std=c++11
    FLAGS = -Wall -O2 -fno-strict-aliasing -ggdb3 -shared
    OUTDIR = lib
    TGT = -DWINDOWS 
else
    CC = g++ -std=c++11
    FLAGS = -Wall -O2 -fno-strict-aliasing -ggdb3
    TGT = -DLINUX
endif

UNIT = -DUNIT_TEST_OWNER
all: unit_geometry unit_lists unit_storage unit_image unit_blob unit_model unit_library

unit_geometry: makefile MWT_Geometry.h MWT_Geometry.cc
	$(CC) $(FLAGS) $(UNIT) $(TGT) -o unit_geometry MWT_Geometry.cc

MWT_Geometry.o: makefile MWT_Geometry.h MWT_Geometry.cc
	$(CC) $(FLAGS) $(TGT) -o MWT_Geometry.o MWT_Geometry.cc
	
unit_lists: makefile MWT_Lists.h MWT_Lists.cc MWT_Geometry.h
	$(CC) $(FLAGS) $(UNIT) $(TGT) -o unit_lists MWT_Lists.cc

MWT_Lists.o: makefile MWT_Lists.h MWT_Lists.cc MWT_Geometry.h
	$(CC) $(FLAGS) $(TGT) -o MWT_Lists.o MWT_Lists.cc
	
unit_storage: makefile MWT_Storage.h MWT_Storage.cc MWT_Lists.h
	$(CC) $(FLAGS) $(UNIT) $(TGT) -o unit_storage MWT_Storage.cc

MWT_Storage.o: makefile MWT_Storage.h MWT_Storage.cc MWT_Lists.h
	$(CC) $(FLAGS) $(TGT) -o MWT_Storage.o MWT_Storage.cc
	
MWT_Image.o: makefile MWT_Storage.h MWT_Geometry.h MWT_Lists.h MWT_Image.h MWT_Image.cc
	$(CC) $(FLAGS) $(TGT) -c -o MWT_Image.o MWT_Image.cc

unit_image: makefile MWT_Storage.h MWT_Geometry.h MWT_Lists.h MWT_Image.h MWT_Image.cc
	$(CC) $(FLAGS) $(UNIT) $(TGT) -o unit_image MWT_Image.cc

MWT_Blob.o: makefile MWT_Storage.h MWT_Geometry.h MWT_Lists.h MWT_Image.h MWT_Blob.h MWT_Blob.cc
	$(CC) $(FLAGS) $(TGT) -c -o MWT_Blob.o MWT_Blob.cc

unit_blob: makefile MWT_Storage.h MWT_Geometry.h MWT_Lists.h MWT_Image.h MWT_Image.o MWT_Blob.h MWT_Blob.cc
	$(CC) $(FLAGS) $(UNIT) $(TGT) -o unit_blob MWT_Blob.cc MWT_Image.o
	
MWT_Model.o: makefile MWT_Storage.h MWT_Geometry.h MWT_Lists.h MWT_Image.h MWT_Image.o MWT_Model.h MWT_Model.cc
	$(CC) $(FLAGS) $(TGT) -c -o MWT_Model.o MWT_Model.cc
	
unit_model: makefile MWT_Storage.h MWT_Geometry.h MWT_Lists.h MWT_Image.h MWT_Image.o MWT_Model.h MWT_Model.cc
	$(CC) $(FLAGS) $(UNIT) $(TGT) -o unit_model MWT_Model.cc MWT_Image.o

unit_library: makefile MWT_Storage.h MWT_Geometry.h MWT_Lists.h MWT_Image.h MWT_Image.o MWT_Blob.h MWT_Blob.o MWT_Model.h MWT_Model.o MWT_Library.h MWT_Library.cc
	$(CC) $(FLAGS) $(UNIT) $(TGT) -o unit_library MWT_Library.cc MWT_Image.o MWT_Blob.o MWT_Model.o

MWT_Library.o: makefile MWT_Storage.h MWT_Geometry.h MWT_Lists.h MWT_Image.h MWT_Blob.h MWT_Model.h MWT_Library.h MWT_Library.cc
	$(CC) $(FLAGS) $(TGT) -c -o MWT_Library.o MWT_Library.cc 

ifeq ($(OS),Windows_NT)
    DLL: makefile MWT_Storage.h MWT_Geometry.h MWT_Lists.h MWT_Image.h MWT_Image.o MWT_Blob.h MWT_Blob.o MWT_Model.h MWT_Model.o MWT_Library.h MWT_Library.o MWT_DLL.h MWT_DLL.cc
	$(CC) -static -static-libstdc++ $(FLAGS) $(TGT) -o $(OUTDIR)/MWT.dll MWT_DLL.cc MWT_Image.o MWT_Blob.o MWT_Model.o MWT_Library.o
endif

test: all
	./unit_geometry
	./unit_lists
	./unit_storage
	./unit_image
	./unit_blob
	./unit_model
	./unit_library

grind: unit_library
	valgrind --leak-check=full --error-exitcode=2 ./unit_library -quiet

clean:
	rm -f unit_geometry unit_lists unit_storage unit_image unit_blob unit_model unit_library
	rm -f test_image.tiff performance_imprint.tiff worm_imprint.tiff worm_noisy.tiff
	rm -f test_blob.log
	rm -f *.o
	rm -f test*.tiff
	rm -f 20071226_105033/*
	rm -f 20071212_130514/*
	rm -fd 20071226_105033 20071212_130514
