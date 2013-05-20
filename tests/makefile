SRCS=$(wildcard *.cpp)
HEADERS=$(SRCS:.cpp=.h)
OBJS=$(SRCS:.cpp=.o)
TARGET=utest

CXXFLAGS := $(shell pkg-config --cflags /home/jninin/Documents/WORK/PROJET/IBEX/OUT3/share/pkgconfig/ibex.pc) -I/home/jninin/Logiciel/cpptest/include/
LIBS	 := -L/home/jninin/Logiciel/cpptest/lib -static $(shell pkg-config --libs  /home/jninin/Documents/WORK/PROJET/IBEX/OUT3/share/pkgconfig/ibex.pc) -lcpptest
 

ifeq ($(DEBUG), yes)
CXXFLAGS := $(CXXFLAGS) -O0 -g -pg -Wall -frounding-math  
else
CXXFLAGS := $(CXXFLAGS) -O3 -DNDEBUG -Wno-deprecated -frounding-math 
endif

$(TARGET):	$(OBJS) 
	$(CXX) -o $(TARGET) $(CXXFLAGS) $(OBJS) $(LIBS)


utest.o : $(SRCS) $(HEADERS) 

%.o :	%.cpp
	$(CXX) $(CPPFLAGS) $(LDFLAGS) $(CXXFLAGS) -c -o $@ $< $(LIBS)

clean:
	rm -f $(OBJS) $(TARGET)
	

