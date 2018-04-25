# Copyright (C) 2006 Simone Furini <sfurini@deis.unibo.it>
#  
# This file is free software; as a special exception the author gives
# unlimited permission to copy and/or distribute it, with or without 
# modifications, as long as this notice is preserved.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY, to the extent permitted by law; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

noinst_HEADERS = functions.h histogram.h

bin_PROGRAMS = village project 

village_SOURCES = village.cc histogram.cc readwrite.cc functions.cc
village_LDADD = -lm
village_LDFLAGS = 
village_DEPENDENCIES = $(includedir)/common.h histogram.h functions.h

project_SOURCES = project.cc functions.cc histogram.cc
project_LDADD = -lm
project_LDFLAGS = 
project_DEPENDENCIES = $(includedir)/common.h histogram.h functions.h
