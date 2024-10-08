//
// Created by fritsche on 07/06/22.
//

#include "data_structures.h"

//ds::Mutations::Mutations() {
//    memset(mutations, 0xFF, m_size);
//}
//
//const size_t ds::Mutations::Size() {
//    size_t size = 0;
//    while (mutations[size] != 0xFF && size < sizeof(mutations)) {
//        size++;
//    }
//    return size;
//}
//
//void ds::Mutations::Add(size_t mutation) {
//    for (size_t i = 0; i < sizeof(mutations); i++) {
//        if (mutations[i] == 0xFF) {
//            mutations[i] = mutation;
//            break;
//        }
//    }
//}

ds::Mutations::Mutations() {
    memset(mutations, 0xFF, size_field);
    mutations[size_field] = 0;
}
void ds::Mutations::Add(size_t mutation) {
    assert(mutations[size_field] < size_field);
    if (Size() == max_size) {
        std::cerr << "Won't add another mutation to the pattern. Size is maxed." << std::endl;
        return;
    }
    mutations[mutations[size_field]++] = mutation;
}

size_t ds::Mutations::Capacity() {
    return sizeof(mutations);
}

bool ds::Mutations::Empty() {
    return mutations[0] == 0xFF;
}

std::string ds::Mutations::ToString(int pos) {
    std::string result = "[";
//    for (auto i = 0; i < sizeof(mutations) && mutations[i] != 0xFF; i++) {
    for (auto i = 0; i < Size(); i++) {
        uint32_t mutation = mutations[i];
        if (i == 0) result += std::to_string(mutation + pos);
        else {
            result += ",";
            result += std::to_string(mutation + pos);
        }
    }
    result += "]";
    return result;
}





ds::SNP::SNP() {}
ds::SNP::SNP(int read_pos, int gene_pos) : read_pos_(read_pos), gene_pos_(gene_pos){};
