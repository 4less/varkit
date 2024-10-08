//
// Created by fritsche on 28/11/2021.
//

#pragma once

#include <sys/types.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <unistd.h>
#include <FFQReader.h>

namespace FFReader {

    static inline int64_t FindNext(char* file, size_t start, size_t end) {
        for (auto i = start; i < end; i++) {
            if (file[i] == '\n') {
                return i;
            }
        }
        return -1;
    }

    static inline int64_t FindNext(char* file, size_t start, size_t end, size_t search_start) {
        //printf("first %c\n", file[search_start]);
        if (file[search_start] == '\n')
            return search_start;

        for (auto osc = start; search_start + osc < end && search_start - osc > start; osc++) {
            if (file[search_start + osc] == '\n')
                return search_start + osc;
            if (file[search_start - osc] == '\n')
                return search_start - osc;
        }

        return FindNext(file, start, end);
    }


    static inline bool ReadRecord(char* file, size_t &start, size_t end, FastxRecord &record) {
        bool ultrafast = false;

        if (start == end || file[start] == '\n')
            return false;

        if (file[start] != '@') {
            std::cerr << "error" << std::endl;
            exit(9);
        }

        static size_t header_len = 0;
        static size_t seq_len = 0;
        static size_t plus_len = 0;
        static size_t qual_len = 0;

        if (ultrafast || header_len == 0) {

            auto header_linebreak = FindNext(file, start, end);
            auto sequence_linebreak = FindNext(file, header_linebreak + 1, end);
            auto plus_linebreak = FindNext(file, sequence_linebreak + 1, end);

            header_len = header_linebreak - start - 1;
            seq_len = sequence_linebreak - header_linebreak - 1;
            plus_len = plus_linebreak - sequence_linebreak - 1;
//            qual_len = qual_start - plus_linebreak - 1;

            record.header.assign(file + start, header_len);
            record.sequence.assign(file + header_linebreak + 1, seq_len);
            record.quality.assign(file + sequence_linebreak + plus_len + 2, seq_len);

            start = plus_linebreak + 1 + seq_len + 1;
        }
        else {
            auto header_linebreak = FindNext(file, start, end, start + header_len + 1);
            auto sequence_linebreak = FindNext(file, header_linebreak + 1, end, header_linebreak + seq_len + 1);
            if (file[sequence_linebreak + 1] != '+')
                sequence_linebreak = FindNext(file, header_linebreak + 1, end);

            if (file[sequence_linebreak + 1] != '+') {
                std::cerr << "error" << std::endl;
                exit(9);
            }
            auto plus_linebreak = FindNext(file, sequence_linebreak + 1, end, sequence_linebreak + plus_len + 1);
//            auto qual_linebreak = FindNext(file, plus_linebreak + 1, end, plus_linebreak + qual_len + 1);

            header_len = header_linebreak - start;
            seq_len = sequence_linebreak - header_linebreak - 1;
            plus_len = plus_linebreak - sequence_linebreak - 1;
//            qual_len = qual_linebreak - plus_linebreak - 1;

            record.header.assign(file + start, header_len);
            record.sequence.assign(file + header_linebreak + 1, seq_len);
            record.quality.assign(file + sequence_linebreak + plus_len + 2, seq_len);
            start = plus_linebreak + 1 + seq_len + 1;
        }

        return true;
    }

    static void ReadFFQ(const char *file_path) {

//        int fd = open(file_path, O_RDONLY, S_IRUSR | S_IWUSR);
//        struct stat sb;
//
//        if (fstat(fd, &sb) == -1) {
//            perror("couldn't get file size. \n");
//        }
//
        FastxRecord record;
        record.format = FORMAT_FASTQ;
//
//        char* file_in_memory = static_cast<char *>(mmap(NULL, sb.st_size, PROT_READ, MAP_PRIVATE, fd, 0));
//
        MapFile *map = new MapFile(file_path);

        size_t start = 0;
        size_t dummy_num = 0;

        while (ReadRecord(map->file_in_memory, start, map->sb.st_size, record)) {
//            std::cout << record.to_string();
            dummy_num += record.sequence.length();
        }
        std::cout << dummy_num << std::endl;
    }

    static void ReadFFQ3(const char *file_path) {


        MapFile *map = new MapFile(file_path);

        FastqRecord record;
        FFQReader reader;

        const size_t block_size = (1024llu * 1024llu * 128);

        size_t dummy_num = 0;

        reader.Set(map);

        while (true) {
            bool ok = false;

            ok = reader.LoadBlock(block_size);
            if (!ok) break;

            while (reader.NextRecord(record)) {
                dummy_num += record.sequence.length();
            }

        }
        map->Close();
        std::cout << dummy_num << std::endl;
    }

    static void ReadFFQ5(const char *file_path) {


        MapFile *map = new MapFile(file_path);


        const size_t block_size = (1024llu * 1024llu * 16llu);
//        const size_t block_size = (1024llu * 1024llu);

        omp_set_num_threads(3);

        size_t dummy_char = 0;
        size_t dummy_num = 0;

#pragma omp parallel
        {
            FastqRecord record;
            FFQReader reader;

#pragma omp critical(reader)
            reader.Set(map);

            while (true) {
                bool ok = false;

#pragma omp critical(reader)
                ok = reader.LoadBlock(block_size);

                if (!ok) break;

                while (reader.NextRecord(record)) {
                    for (auto i = 0; i < record.sequence.length(); i++) {
#pragma omp atomic
                        dummy_char += record.sequence[i];
                    }
#pragma omp atomic
                    dummy_num += record.sequence.length();
                }
            }
        }
        map->Close();
        std::cout << dummy_num << std::endl;
        std::cout << dummy_char << std::endl;
    }

    static void ReadFFQ4(const char *file_path) {
        const size_t block_size = (1024 * 1024);
        ifstream is(file_path, ios::in);

        size_t dummy_char = 0;
        size_t dummy_num = 0;

        omp_set_num_threads(3);

#pragma omp parallel
        {
            BufferedFastxReader reader;
            FastxRecord record;

            while (true) {
                bool ok = false;

#pragma omp critical(reader)
                ok = reader.LoadBlock(is, block_size);

                if (!ok) break;

                // Read records from datablock
                while (true) {
                    auto valid_fragment = reader.NextSequence(record);
                    if (!valid_fragment) break;

                    for (auto i = 0; i < record.sequence.length(); i++) {
#pragma omp atomic
                        dummy_char += record.sequence[i];
                    }
#pragma omp atomic
                    dummy_num += record.sequence.length();
                } // End Read Block
            } // End Reader
        } // Close OMP parallel
        is.close();

        std::cout << dummy_num << std::endl;
        std::cout << dummy_char << std::endl;
    }

    static void ReadFFQ2(const char *file_path) {
        ifstream is(file_path, ios::in);

        const size_t block_size = (1024 * 1024);

        BufferedFastxReader reader;
        FastxRecord record;

        size_t dummy_num = 0;

        while (true) {
            bool ok = false;

            ok = reader.LoadBlock(is, block_size);
            if (!ok) break;

            // Read records from datablock
            while (true) {
                auto valid_fragment = reader.NextSequence(record);
                if (!valid_fragment) break;

                dummy_num += record.sequence.length();
//                std::cout << record.to_string();
            }
        }
        is.close();
        std::cout << dummy_num << std::endl;
    }

    static void test() {
        std::string file = "/usr/users/QIB_fr017/fritsche/Projects/data/camisim/short_read/2017.12.04_18.45.54_sample_1/reads/anonymous_reads.fq";
//        std::string file = "/usr/users/QIB_fr017/fritsche/Projects/data/reads/dummy.fq";
//        std::string file = "/usr/users/QIB_fr017/fritsche/Projects/data/camisim/short_read/2017.12.04_18.45.54_sample_1/reads/anon_sub.fq";

//
//        std::cout << "New Reader" << std::endl;
//        Benchmark bm3("New Reader");
//        bm3.start();
//        ReadFFQ3(file.c_str());
//        bm3.stop();



        std::cout << "Old Reader MT" << std::endl;
        ds::Benchmark bm4("Old Reader MT");
        bm4.start();
        ReadFFQ4(file.c_str());
        bm4.stop();

        std::cout << "New Reader MMT" << std::endl;
        ds::Benchmark bm5("New Reader MMT");
        bm5.start();
        ReadFFQ5(file.c_str());
        bm5.stop();

//        std::cout << "Old Reader" << std::endl;
//        Benchmark bm2("Old Reader");
//        bm2.start();
//        ReadFFQ2(file.c_str());
//        bm2.stop();

//        std::cout << "New Method" << std::endl;
//        Benchmark bm1("New Method");
//        bm1.start();
//        ReadFFQ(file.c_str());
//        bm1.stop();

//        bm1.printResults();
//        bm2.printResults();
//        bm3.printResults();
        bm4.printResults();
        bm5.printResults();
    }
}