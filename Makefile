bin/cvg : bin/cvg.o bin/common.o bin/ReadMapping.o bin/ReadMappingReader.o bin/ReadSequence.o bin/String.o
	g++ bin/cvg.o bin/common.o bin/ReadMapping.o bin/ReadMappingReader.o bin/ReadSequence.o bin/String.o -o bin/cvg

bin/cvg.o : cvg.cpp
	g++ -c cvg.cpp -o bin/cvg.o

bin/common.o : common/common.cpp
	g++ -c common/common.cpp -o bin/common.o

bin/ReadMapping.o : common/ReadMapping.cpp
	g++ -c common/ReadMapping.cpp -o bin/ReadMapping.o

bin/ReadMappingReader.o : common/ReadMappingReader.cpp
	g++ -c common/ReadMappingReader.cpp -o bin/ReadMappingReader.o

bin/ReadSequence.o : common/ReadSequence.cpp
	g++ -c common/ReadSequence.cpp -o bin/ReadSequence.o

bin/String.o : common/String.cpp
	g++ -c common/String.cpp -o bin/String.o

