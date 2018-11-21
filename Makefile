ESPRESSOPPDIR = /home/niklas/Documents/Studium/BA/ba_mamico/espressopp
.PHONY: espressopp

ESPRESSOPPFLAGS=-DCMAKE_BUILD_

CXX = mpicxx
CFLAGS= -g3 -O0 -Wall --std=c++11

RM = rm -rf

INCLUDES = -I/usr/include/python2.7 \
           -I$(ESPRESSOPPDIR) \
           -I$(ESPRESSOPPDIR)/src \
           -I$(ESPRESSOPPDIR)/src/include \
           -I$(ESPRESSOPPDIR)/contrib/boost \
           -I/home/wittmer/ba/mamico/coupling/interface/impl/Espresso++ \
           -I/home/wittmer/.local/share/include/xdrfile \

EPP_LIB_PATH = -L$(ESPRESSOPPDIR) -L/home/wittmer/.local/share/lib
EPP_LIB = :_espressopp.so

LDFLAGS = -lpython2.7 -l$(EPP_LIB) 

OBJ = ../src/io/DumpXTC.o \
      ../../coupling/interface/impl/Espresso++/ReflectingBC.o

LJSimple: LennardJonesSimple.o $(OBJ)
	$(CXX) $(CFLAGS) $^ $(INCLUDES) $(EPP_LIB_PATH) $(LDFLAGS) -o $@ 

LJAdress: LennardJonesAdress.o LJ_include.hpp
	$(CXX) $(CFLAGS) $< $(INCLUDES) $(EPP_LIB_PATH) $(LDFLAGS) -o $@

%.o: %.cpp
	$(CXX) $(CFLAGS) $< $(INCLUDES) $(EPP_LIB_PATH) $(LDFLAGS) -c -o $@

espressopp:
	$(MAKE) -C $(ESPRESSOPPDIR)

clean:
	$(RM) LJSimple
	$(RM) LJAdress
	$(RM) $(OBJ)
	$(RM) LennardJonesSimple.o
	$(RM) LennardJonesAdress.o
