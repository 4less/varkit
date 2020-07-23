cmake_minimum_required(VERSION 3.15)
project(varkit)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS "-pthread")

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")
link_directories(/home/joachim/library_cpp/bhm/lib)

# include src dirs
link_directories(src/io src/kmer_processing src/taxonomy src/variant_predictor src/test lib src/utils)
include_directories(src/io src/test lib/ src/kmer_processing src/utils src/taxonomy)

add_executable(varkit
        src/main.cpp
        src/VarkitExecutor.h
        src/io/MetaDataDB.cpp
        src/io/MetaDataDB.h
        src/io/FastaReader.cpp
        src/io/FastaReader.h
        src/io/FastaReader2.cpp
        src/io/FastaReader2.h
        src/utils/ShapeUtils.h
        src/utils/Utils.h
        src/test/test.h
        src/OptionsContainer.h
        src/kmer_processing/KmerIterator.h
        src/kmer_processing/KmerIterator.cpp
        src/kmer_processing/SpacedKmerIterator.h
        src/kmer_processing/SpacedKmerIterator.cpp
        src/kmer_processing/KmerProcessor.h
        src/utils/KmerUtils.h src/test/iotest.h
        src/taxonomy/Taxonomy.h)


include(FindPkgConfig)

# include BitHashMap
find_package(BHM REQUIRED)
target_include_directories(${PROJECT_NAME} PUBLIC ${BHM_INCLUDE_DIR})
target_link_libraries(${PROJECT_NAME} LINK_PUBLIC ${BHM_LIBRARY})

# include CWCSV
find_package(CWCSV REQUIRED)
target_include_directories(${PROJECT_NAME} PUBLIC ${CWCSV_INCLUDE_DIR})

include_directories(/home/joachim/library_cpp/robin/include)



