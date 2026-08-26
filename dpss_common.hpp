#pragma once
#include <algorithm>
#include <atomic>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <new>
#include <expected>
#include <format>
#include <fstream>
#include <functional>
#include <iostream>
#include <queue>
#include <random>
#include <ranges>
#include <sstream>
#include <string>
#include <vector>

using std::size_t;
using std::vector;

namespace allocstat {
    static std::atomic<unsigned long long> count{0};
    static std::atomic<unsigned long long> bytes{0};
    static inline unsigned long long allocs() { return count.load(std::memory_order_relaxed); }
    static inline unsigned long long nbytes() { return bytes.load(std::memory_order_relaxed); }
}
void* operator new(std::size_t sz) {
    allocstat::count.fetch_add(1, std::memory_order_relaxed);
    allocstat::bytes.fetch_add(sz, std::memory_order_relaxed);
    if (void* p = std::malloc(sz ? sz : 1)) return p;
    throw std::bad_alloc();
}
void* operator new[](std::size_t sz) { return ::operator new(sz); }
void operator delete(void* p) noexcept { std::free(p); }
void operator delete(void* p, std::size_t) noexcept { std::free(p); }
void operator delete[](void* p) noexcept { std::free(p); }
void operator delete[](void* p, std::size_t) noexcept { std::free(p); }

namespace deadline {
    using Clock = std::chrono::steady_clock;
    static Clock::time_point limit;
    static bool armed = false;
    static bool timedOut = false;
    static inline void arm(double seconds) {
        timedOut = false;
        if (seconds > 0.0) {
            limit = Clock::now() +
                    std::chrono::duration_cast<Clock::duration>(
                        std::chrono::duration<double>(seconds));
            armed = true;
        } else {
            armed = false;
        }
    }
    static inline bool hit() {
        if (armed && Clock::now() >= limit) { timedOut = true; return true; }
        return false;
    }
}
