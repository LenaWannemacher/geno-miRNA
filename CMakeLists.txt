cmake_minimum_required(VERSION 3.14)
project(mirna)

set(CMAKE_CXX_STANDARD 17)

add_library(SHARED
        json.hpp
        InvalidCharacter.cpp
        InvalidCharacter.h)

add_executable(mirna
        main.cpp
        InvalidCharacter.cpp
        InvalidCharacter.h
        SeqAlignment.h
        algorithm.h
        output.h
        intervallTree.h
        Candidate.h
        GenomeProcessor.h
        GenomeProcessor.cpp
        GFFReader.h
        GFFReader.cpp
        Matrix.h
        StringMatrix.h
        miRNA.h
        miRNA.cpp)