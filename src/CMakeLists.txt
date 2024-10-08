project(varkit_lib)
set(CMAKE_CXX_STANDARD 20)



set(CMAKE_EXE_LINKER_FLAGS "-lz")
set(CMAKE_CXX_FLAGS_DEBUG "-pthread -ggdb -g -std=c++2a -fopenmp -Wall")
set(CMAKE_CXX_FLAGS_RELEASE "-std=c++2a -DNDEBUG -pthread -O3 -fopenmp -Werror=unknown-pragmas") #-march=native


# This creates the target and must be specified before using target_include etc
add_library(varkit_lib
        modes/build.h
        modes/classify.h
        modes/profile.h

        compact_map/multimap.h
        minimap/minimap.h
        OptionsContainer.h

#        io/MetaDataDB.cpp
#        io/MetaDataDB.h
#        io/FastaReader.cpp
#        io/FastaReader.h
#        io/FastxReader.cpp
#        io/FastxReader.h
#        io/FastaReader2.cpp
#        io/FastaReader2.h
#        io/BufferedOutput.h
#        io/FFReader.h
#        io/FFQReader.h
#        io/VCFHandler.h
#        io/GenomeLoader.h
#        io/BufferedOutput.cpp
#
#
#
#        kmer_processing/KmerIterator.h
#        kmer_processing/KmerIterator.cpp
#        kmer_processing/SpacedKmerIterator.h
#        kmer_processing/SpacedKmerIterator.cpp
#        kmer_processing/KmerProcessor.h
#        kmer_processing/SimpleKmerPutter.h
#        kmer_processing/SimpleKmerIterator.h
#
#        compact_map/compact_map.h
#        compact_map/compact_map.cpp
#
#        variant_predictor/ShapeBuilder.h
#        variant_predictor/ShapeBuilder.cpp
#        variant_predictor/CombinationIterator.h
#        variant_predictor/ShapeHistory.h
#        variant_predictor/ShapeHistory.cpp
#        variant_predictor/SNPPredictor.h
#        variant_predictor/SNPPredictor.cpp
#        variant_predictor/IOHandler.h
#
        utils/thread_pool.h
#        utils/Progress.h
#        utils/Progress.cpp
#        utils/BinaryClassifierEvaluator.h
#        utils/tiny_tree.h
#        utils/version.h
#        utils/KmerUtils.h
#        utils/KmerBuffer.h
#        utils/ShapeUtils.h
#        utils/Utils.h
#        utils/progress_bar.h
#
#        taxonomy/TaxonomyNew.cpp
#        taxonomy/TaxonomyNew.h
#        taxonomy/Test.h
#        taxonomy/CustomTaxonomy.h
#        taxonomy/Taxonomy.h
#        taxonomy/TestTaxonomy.h
        taxonomy/ReadClassifier.h

#        classify/classify_scripts.h
#
#        test/VCFTest.h
#        test/kmer_iterator_test.h
#        test/test2.h
#        test/test.h
#        test/iotest.h
#
#
#        classify/AbundanceEstimator.h
#        classify/classification.h
#        classify/classification.cpp
#        classify/AbundanceEstimator.cpp
#        classify/strain_handler.h
#        classify/strain_handler.cpp
#
#
#        data_structures/data_structures.h
#        data_structures/data_structures.cpp
#        varkit_config.h
#        constants.h


        #------
#        align/Aligner.h
#        align/Aligner.cpp
        #        utils/ThreadPool.h
        #        utils/semaphore.h
        #        utils/queue.h
        ##        build/build_scripts.h
        ##        taxonomy/HitsClassifier.cpp
        ##        taxonomy/HitsClassifier.h
        ##        variant_predictor/VariantInferer.h
        ##        variant_predictor/VariantInferer.cpp
        ##        VarkitExecutor.h
)


function(get_all_targets var)
    set(targets)
    get_all_targets_recursive(targets ${CMAKE_CURRENT_SOURCE_DIR})
    set(${var} ${targets} PARENT_SCOPE)
endfunction()

macro(get_all_targets_recursive targets dir)
    get_property(subdirectories DIRECTORY ${dir} PROPERTY SUBDIRECTORIES)
    foreach(subdir ${subdirectories})
        get_all_targets_recursive(${targets} ${subdir})
    endforeach()

    get_property(current_targets DIRECTORY ${dir} PROPERTY BUILDSYSTEM_TARGETS)
    list(APPEND ${targets} ${current_targets})
endmacro()

get_all_targets(all_targets)
message("All targets: ${all_targets}")


###########################################################################################
# Now target varkit_lib exists
###########################################################################################

# link_directories makes files in these directories includable without explicitly stating the dir
# include "test/kmer_iterator_test.h" -> include "kmer_iterator_test.h"
#link_directories(varkit_lib
#        PUBLIC io/
#        kmer_processing/
#        taxonomy/
#        variant_predictor/
#        test/
#        lib/
#        utils/)

target_link_directories(varkit_lib
        PUBLIC
        io/
        kmer_processing/
        taxonomy/
        variant_predictor/
        test/
        lib/
        utils/)

target_include_directories(varkit_lib PUBLIC
        ../lib/
        build/
        classify/
        compact_map/
        io/
        kmer_processing/
        minimap/
        modes/
        taxonomy/
        test/
        utils/
        variant_predictor/
        data_structures)


# Adding this to the path makes find_package find the required BHM library
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/../cmake")
#link_directories(/usr/users/QIB_fr017/fritsche/library_cpp/bhm/lib)


#message(STATUS "BHM include dir: ${BHM_INCLUDE_DIR}")
message(STATUS "BHM library: ${BHM_LIBRARY}")
message(STATUS "CSV library: ${CWCSV_INCLUDE_DIR}")

# Directly include folders to make them "includable" with '#include "robin_map.h"'
# require BitHashMap
#find_package(BHM REQUIRED)
find_package(ZLIB REQUIRED)

# include CWCSV
find_package(CWCSV REQUIRED)


# Directly include folders to make them "includable" with '#include "robin_map.h"'
# generic include_directories should link with all targets in "dir"
# however this did not work, so I had to change to target_include_directories

# Both needed (for whatever reason)
#target_include_directories(varkit_lib PUBLIC ${BHM_INCLUDE_DIR})
target_link_libraries(varkit_lib LINK_PUBLIC ${BHM_LIBRARY})

target_include_directories(varkit_lib PUBLIC /usr/users/QIB_fr017/fritsche/library_cpp/robin)
target_include_directories(varkit_lib PUBLIC /usr/users/QIB_fr017/fritsche/library_cpp/tsl)

# Necessary to link to stdc++fs (filesystem experimental support)
target_link_libraries(varkit_lib LINK_PUBLIC ZLIB::ZLIB stdc++fs)


target_include_directories(varkit_lib PUBLIC ${CWCSV_INCLUDE_DIR})
#target_include_directories(${PROJECT_NAME} PUBLIC ${CWCSV_INCLUDE_DIR})
set_target_properties(varkit_lib PROPERTIES LINKER_LANGUAGE CXX)