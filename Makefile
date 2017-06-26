# Minimalist makefile

CXXFLAGS=-std=c++14 -pthread -Wno-unused-function -DNDEBUG -O2 -Wall -IExternalSoftware -ISrc
LDFLAGS=-lm -llapacke -llapack -lcblas -lblas 
SRC= $(wildcard Src/*.cpp)
OBJ= $(SRC:.cpp=.o)

all: jointDeconvolution generateSynthetic

libjointDeconv.a: $(OBJ)
	ar rcs libjointDeconv.a $^

jointDeconvolution: jointDeconvolution.cpp libjointDeconv.a
	$(CXX) $(CXXFLAGS) -o $@ $< libjointDeconv.a $(LDFLAGS)

generateSynthetic: generateSynthetic.cpp libjointDeconv.a
	$(CXX) $(CXXFLAGS) -o $@ $< libjointDeconv.a $(LDFLAGS)

clean:
	rm $(OBJ) jointDeconvolution generateSynthetic libjointDeconv.a
