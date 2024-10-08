//
// Created by fritsche on 28/11/2021.
//

#ifndef VARKIT_FFQREADER_H
#define VARKIT_FFQREADER_H

#include <sys/types.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <unistd.h>
#include <iostream>
#include <bits/shared_ptr.h>
#include "FastxReader.h"

struct FastqRecord {
    FileFormat format;
    std::string header = "";  // header line, including @/>, but not newline
    std::string id = "";      // from first char. after @/> up to first whitespace
    std::string sequence = "";
    std::string quality = "";   // only meaningful for FASTQ seqs

    std::string &ToString() {
        record_str.assign(header);
        record_str.append("\n");
        record_str.append(sequence);
        record_str.append("\n+\n");
        record_str.append(quality);
        record_str.append("\n");
        return record_str;
    }

private:
    std::string record_str;
};

class MapFile {
public:
    char* file_in_memory;
    size_t current_pos = 0;
    int fd;
    struct stat sb;

    MapFile(const char* file_path) {
        fd = open(file_path, O_RDONLY, S_IRUSR | S_IWUSR);
        if (fstat(fd, &sb) == -1) {
            perror("couldn't get file size. \n");
        }
        file_in_memory = (char*) mmap(NULL, sb.st_size, PROT_READ, MAP_PRIVATE, fd, 0);
    }

    void Close() {
        close(fd);
        munmap(file_in_memory, Size());
    }

    size_t Size() {
        return sb.st_size;
    }
};

class FFQReader {
private:
    size_t block_start_;
    size_t block_end_ = 0;

    MapFile *map;
    char* file;

    static inline int64_t FindNext(char* file, size_t start, size_t end) {
        for (auto i = start; i < end; i++) {
            if (file[i] == '\n')
                return i;
        }
        std::cout << "return -1" << std::endl;
        exit(15);
        return -1;
    }

public:
    size_t current_position_ = 0;

    void Set(MapFile *map) {
        this->map = map;
        this->file = map->file_in_memory;
        this->block_start_ = 0;
        this->current_position_ = map->current_pos;
    }

    void Print(size_t start, size_t len) {
        std::cout << "print start: " << start << std::endl;
        for (size_t i = start; i < start+len; i++) {
            std::cout << file[i];
        }
        std::cout << std::endl;
        std::cout << "___" << std::endl;
    }

    bool LoadBlock(const size_t target_block_size) {

//        if (file[current_position_] != '@' ||
//            (current_position_ > 0 && file[current_position_ - 1] != '\n') ||
//            (file[current_position_+1] != 'S')) {
//            std::cout << omp_get_thread_num() << " problem at: " << map->current_pos << std::endl;
//            std::cout << omp_get_thread_num() << " problem at reader curpos: " << current_position_ << std::endl;
//            std::cout << omp_get_thread_num() << " problem at reader block end: " << block_end_ << std::endl;
//            for (auto i = current_position_; i < current_position_ + 100; i++) {
//                std::cout << file[i];
//            }
//            std::cout << "_________" << std::endl;
//
//            exit(48);
//        }




        if (map->current_pos >= map->Size()) {
            return false;
        } else {

//            if (file[map->current_pos] != '@' ||
//                (map->current_pos > 0 && file[map->current_pos - 1] != '\n') ||
//                (file[map->current_pos+1] != 'S')) {
//                std::cout << omp_get_thread_num() << " problem at reader curpos: " << current_position_ << std::endl;
//                std::cout << omp_get_thread_num() << " problem at reader block end: " << block_end_ << std::endl;
//
//                std::cout << omp_get_thread_num() << " problem at: " << map->current_pos << std::endl;
//                for (auto i = map->current_pos; i < map->current_pos + 100; i++) {
//                    std::cout << file[i];
//                }
//                std::cout << "_________" << std::endl;
//
//                exit(49);
//            }

            // Move current pos to global pos in file
            this->current_position_ = map->current_pos;
            if (map->current_pos + target_block_size >= map->Size()) {
                // Block end needs to be complete file end
                this->block_end_ = map->Size();
            } else {
                // Block end is current_pos + block_size

                this->block_end_ = this->current_position_ + target_block_size;

                auto end_save = block_end_;

                // Find header start
                while (true) {
                    if (file[this->block_end_] != '@' || file[this->block_end_ - 1] != '\n') {
                        this->block_end_--;
                        continue;
                    }
                    auto j = this->block_end_ - 2;
                    while (file[j] != '\n') j--;
                    if (file[j+1] == '+') {
                        this->block_end_ = j - 1;
                        continue;
                    } else {
                        break;
                    }
                }
//                std::cout << "print block end" << std::endl;
//                Print(block_end_, 500);


//                if (file[this->block_end_] != '@' || file[this->block_end_-1] != '\n') {
//                    exit(34);
//                }
            }
            // move on global pos to end of block
            map->current_pos = block_end_;
        }

//        std::cout << omp_get_thread_num() << " LoadBlock " << current_position_ << " - " << block_end_ << std::endl;

        return true;
    }

    bool NextRecord(FastqRecord &record) {

        if (current_position_ == block_end_) {
            return false;
        }
//        if (file[current_position_] == '\n') {
//            exit(99);
//        }

//        if (file[current_position_] != '@' ||
//                (current_position_ > 0 && file[current_position_ - 1] != '\n') ||
//                (file[current_position_+1] != 'S')) {
//            std::cout << omp_get_thread_num() << " problem at: " << current_position_ << std::endl;
//            for (auto i = current_position_; i < current_position_ + 100; i++) {
//                std::cout << file[i];
//            }
//            std::cout << "_" << std::endl;
//
//            std::cout << "previous record: " << std::endl;
//            std::cout << "length: " << record.sequence.length() << std::endl;
//            std::cout << record.sequence << std::endl;
//            std::cout << record.quality << std::endl;
//            std::cout << "___" << std::endl;
//            std::cout << record.ToString();
//            exit(54);
//        }

        auto header_linebreak = FindNext(file, current_position_, block_end_);
        auto sequence_linebreak = FindNext(file, header_linebreak + 1, block_end_);
        auto plus_linebreak = FindNext(file, sequence_linebreak + 1, block_end_);

//        if (file[sequence_linebreak+1] != '+') {
//            exit(56);
//        }
//        if (file[current_position_] != '@') {
//            exit(57);
//        }


        auto header_len = header_linebreak - current_position_ - 1;
        auto seq_len = sequence_linebreak - header_linebreak - 1;
        auto plus_len = plus_linebreak - sequence_linebreak - 1;

//        if (file[plus_linebreak + seq_len + 1] != '\n') {
//            exit(58);
//        }
//        if (file[plus_linebreak + seq_len + 2] != '@') {
//            exit(59);
//        }

//        if (plus_linebreak + seq_len + 2 == 402652612) {
//            std::cout << "awaiting break..." << std::endl;
//        }

//        if (plus_linebreak + seq_len + 1 > block_end_) {
//            return false;
//        }


        record.header.assign(file + current_position_, header_len);
        record.sequence.assign(file + header_linebreak + 1, seq_len);
        record.quality.assign(file + sequence_linebreak + plus_len + 2, seq_len);

        current_position_ = plus_linebreak + seq_len + 2;

        if (current_position_ >= map->Size()) {
            return true;
        }

//        if (file[current_position_] != '@' ||
//            (current_position_ > 0 && file[current_position_ - 1] != '\n') ||
//            (file[current_position_+1] != 'S') ||
//            record.sequence.size() < 10) {
//            std::cout << omp_get_thread_num() << " problem at: " << current_position_ << std::endl;
//            std::cout << (file[current_position_] != '@') << std::endl;
//            std::cout << (current_position_ > 0 && file[current_position_ - 1] != '\n') << std::endl;
//            std::cout << (file[current_position_+1] != 'S') << std::endl;
//            std::cout << (record.sequence.size() < 10) << std::endl;
//            std::cout << "___________100 chars in file:" << std::endl;
//            for (size_t i = current_position_; i < current_position_ + 100; i++) {
//                std::cout << file[i];
//            }
//            std::cout << "____end____" << std::endl;
//
//            std::cout << "previous record: " << std::endl;
//            std::cout << "length: " << record.sequence.length() << std::endl;
//            std::cout << record.sequence << "." << std::endl;
//            std::cout << record.quality << "." << std::endl;
//            std::cout << "___" << std::endl;
//            std::cout << record.ToString();
//            exit(55);
//        }
        return true;
    }

};


#endif //VARKIT_FFQREADER_H
