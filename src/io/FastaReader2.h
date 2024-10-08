//
// Created by joachim on 06/07/2020.
//

#ifndef VARKIT_FASTAREADER2_H
#define VARKIT_FASTAREADER2_H
#define BUF_SIZE 1 << 17

#include <fstream>

using namespace std;

namespace FastaReader2 {

    static int pos = BUF_SIZE;
    static char buf[BUF_SIZE];
    
    inline char nextch(FILE *fin) {
        if (pos == BUF_SIZE) fread(buf, BUF_SIZE, 1, fin), pos = 0;
        return buf[pos++];
    }
    
    inline int read(FILE *fin) {
        thread_local static int pos = 0;
        thread_local static char ch;
        while ((ch = nextch(fin)) != '>');
        int x = ch - '0';
        while (isdigit(ch = nextch(fin))) x = 10 * x + ch - '0';
        return x;
    }
}


#endif //VARKIT_FASTAREADER2_H
