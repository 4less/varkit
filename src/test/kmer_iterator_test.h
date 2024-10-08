//
// Created by fritsche on 11/03/2021.
//

#pragma once

#include <string>
#include <SimpleKmerIterator.h>
#include "ShapeUtils.h"

namespace KmerIteratorTest {
    static void test() {
        std::string ref = "ACGTACTACTGATCACTATCTAGCTAGTCGTAGC";
        std::string shape_str = "X__XX_X_XX__X";

        std::cout << ref << std::endl;
        std::cout << shape_str << std::endl;

        bool* shape = ShapeUtils::GetShape(shape_str);

        uint64_t key;
        size_t k = 7;

        SimpleKmerIterator iterator(k, shape, shape_str.length());

        FastxRecord record;
        record.sequence = ref;

        iterator.SetRecord(record);

        auto pos = 0;
        while (iterator.HasNext()) {
            iterator.operator()(key);
//            std::cout << ref << std::endl;
//            std::cout << std::string(pos, ' ') << shape_str << std::endl;
//            std::cout << std::string(pos, ' ') << KmerUtils::ToString(key, 2*k) << std::endl;
//            std::cout << "hasnext: " << (iterator.HasNext()) << std::endl;
            ++pos;
        }
    }
}