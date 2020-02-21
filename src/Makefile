EXECUTABLE=StLiter

CXX = g++

CFLAGS = -lpthread

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
endif
ifeq ($(UNAME_S),Darwin)
	CFLAGS += -I/opt/local/include
endif

$(EXECUTABLE): Main.cpp Basic.o Hash.o litool.o Combine_Kmers.o Extend_Kmer.o Final_Extend.o Being_Short.o Builder.o BplusTreeBit.o 
	$(CXX) $^ -o $@ $(CFLAGS) 

Basic.o: Basic.cpp
	$(CXX) -c $< -o $@

Hash.o: Hash.cpp #Basic.h BplusTree.h
	$(CXX) -c $< -o $@

litool.o: litool.cpp #Basic.h BplusTree.h BplusTreeBit.h Final_Extend.h
	$(CXX) -c $< -o $@

BplusTreeBit.o: BplusTreeBit.cpp 
	$(CXX) -c $< -o $@

Combine_Kmers.o: Combine_Kmers.cpp #Basic.h
	$(CXX) -c $< -o $@	

Extend_Kmer.o: Extend_Kmer.cpp #Basic.h BplusTree.h Hash.h Final_Extend.h litool.h Builder.h
	$(CXX) -c $< -o $@	

Final_Extend.o: Final_Extend.cpp #Basic.h BplusTree.h BplusTreeBit.h Hash.h litool.h
	$(CXX) -c $< -o $@	

Being_Short.o: Being_Short.cpp #Basic.h
	$(CXX) -c $< -o $@

Builder.o: Builder.cpp #Basic.h BplusTree.h
	$(CXX) -c $< -o $@

.PHONY : clean

clean:
	rm *.o $(EXECUTABLE)
