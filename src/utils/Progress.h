//
// Created by fritsche on 29/04/2021.
//

#pragma once

#include <cstddef>
#include <mutex>
#include "progress_bar.h"

class ProgressInterface {
public:
    virtual void Reset(size_t progress_size) = 0;
    virtual void Update(size_t new_progress) = 0;
    virtual void UpdateAdd(size_t new_progress) = 0;
};

class ProgressSingleton : public ProgressInterface {
    static ProgressSingleton *instance;
    ProgressBar bar;
    std::mutex lock;

public:
    ProgressSingleton() {};
    ~ProgressSingleton() {
        delete instance;
    }

    static ProgressSingleton *GetInstance() {
        if (!instance)
            instance = new ProgressSingleton();
        return instance;
    }

    void Reset(size_t progress_size) override {
        std::lock_guard<std::mutex> progress_lock(lock);
        bar.reset(progress_size);
    }

    void Update(size_t new_progress) override {
        std::lock_guard<std::mutex> progress_lock(lock);
        bar.Update(new_progress);
    }

    void UpdateAdd(size_t new_progress) override {
        std::lock_guard<std::mutex> progress_lock(lock);
        bar.UpdateAdd(new_progress);
    }
};