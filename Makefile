# Minimalist makefile

#CXXFLAGS=-std=c++14 -pthread -Wno-unused-function -DNDEBUG -O2 -Wall -IExternalSoftware -ISrc -Wno-format-truncation
CXXFLAGS=-std=c++14 -pthread -Wno-unused-function -g -Wall -IExternalSoftware -ISrc -Wno-format-truncation
LDFLAGS=-lm -llapacke -llapack -lblas 
SRC= $(wildcard Src/*.cpp)
OBJ= $(SRC:.cpp=.o)

all: jointDeconvolution generateSynthetic sequencialDeconvolution

libjointDeconv.a: $(OBJ)
	ar rcs libjointDeconv.a $^

jointDeconvolution: jointDeconvolution.cpp libjointDeconv.a
	$(CXX) $(CXXFLAGS) -o $@ $< libjointDeconv.a $(LDFLAGS)

sequencialDeconvolution: sequencialDeconvolution.cpp libjointDeconv.a
	$(CXX) $(CXXFLAGS) -o $@ $< libjointDeconv.a $(LDFLAGS)

generateSynthetic: generateSynthetic.cpp libjointDeconv.a
	$(CXX) $(CXXFLAGS) -o $@ $< libjointDeconv.a $(LDFLAGS)

clean:
	rm $(OBJ) jointDeconvolution generateSynthetic libjointDeconv.a
