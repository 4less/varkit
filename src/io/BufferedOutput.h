//
// Created by fritsche on 28/11/2021.
//

#ifndef VARKIT_BUFFEREDOUTPUT_H
#define VARKIT_BUFFEREDOUTPUT_H

//#include <bits/stdc++.h>
#include <cstddef>
#include <fstream>
#include <assert.h>

using namespace std;

enum BufferedOutputUnit {
    N, MB, GB
};

template<typename T>
class BufferedOutput {
    T* buffer_;

    size_t buffer_fill_ = 0;

    size_t buffer_size_;
    static constexpr size_t entry_size_ = sizeof(T);


public:
    BufferedOutput(size_t buffer_size) {
        buffer_size_ = buffer_size;
        buffer_ = new T[buffer_size_];
    }

    BufferedOutput(size_t buffer_size, BufferedOutputUnit bsu) {
        switch (bsu) {
            case N:
                buffer_size_ = buffer_size;
                break;
            case MB:
                buffer_size_ = (buffer_size*1024*1024) / sizeof(T);
                break;
            case GB:
                buffer_size_ = (buffer_size*1024*1024*1024) / sizeof(T);
                break;
        }
//        std::cout << "Allocate with: " << buffer_size_ << " -> " << (sizeof(T) * buffer_size_)/(1024*1024) << " MB" << std::endl;

        buffer_ = new T[buffer_size_];
    }

    ~BufferedOutput() {
        delete[] buffer_;
    }

    inline bool Write(T &item) {
        assert(buffer_fill_ + entry_size_ < buffer_size_);
        memcpy(buffer_ + buffer_fill_, &item, entry_size_);
        buffer_fill_++;
        return buffer_fill_ < buffer_size_;
    }

    inline void Write(std::ostream &ofs) {
        if (buffer_fill_ != 0) {
//            std::cout << "write snps" << buffer_fill_ << std::endl;
            ofs.write(reinterpret_cast<char *>(buffer_), buffer_fill_ * sizeof(T));
            buffer_fill_ = 0;
        }
    }
};


#endif //VARKIT_BUFFEREDOUTPUT_H
