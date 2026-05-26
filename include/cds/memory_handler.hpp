//
// Created by Diego Diaz on 2/10/20.
//

#ifndef CDS_MEMORY_HANDLER_HPP
#define CDS_MEMORY_HANDLER_HPP

#ifdef __linux__
#include <malloc.h>
#endif

struct mem {

    template<class T>
    static T * allocate(const size_t n) {
        if (n == 0) {return nullptr;}
        if (n > static_cast<size_t>(-1) / sizeof(T)) {
            throw std::bad_array_new_length();
        }
        void * const pv = malloc(n * sizeof(T));
        if (!pv) { throw std::bad_alloc(); }
        return static_cast<T *>(pv);
    }

    template<class T>
    static T * reallocate(T* old_ptr, const size_t n) {
        if (n == 0) {return nullptr;}
        if (n > static_cast<size_t>(-1) / sizeof(T)) {
            throw std::bad_array_new_length();
        }
        void * const pv = realloc(old_ptr, n * sizeof(T));
        if (!pv) { throw std::bad_alloc(); }
        return static_cast<T *>(pv);
    }

    template<class T>
    static void deallocate(T *& p) {
        if (p!= nullptr) {
            free(p);
#ifdef __linux__
            malloc_trim(0);
#endif
            p=nullptr;
        }
    }
};

/* Obtained from
 * https://stackoverflow.com/questions/36517825/is-stephen-lavavejs-mallocator-the-same-in-c11/36521845#36521845
 */
template <class T> struct mallocator {

    typedef T value_type;

    mallocator() noexcept = default; // default ctor not required

    template <class U> explicit mallocator(const mallocator<U>&) noexcept { }

    template <class U>
    bool operator==(const mallocator<U>&) const noexcept {
        return true;
    }

    template <class U> bool operator!=(const mallocator<U>&) const noexcept {
        return false;
    }

    static T * allocate(const size_t n) {
        if (n == 0) {return nullptr;}
        if (n > static_cast<size_t>(-1) / sizeof(T)) {
            throw std::bad_array_new_length();
        }
        void * const pv = malloc(n * sizeof(T));
        if (!pv) { throw std::bad_alloc(); }
        return static_cast<T *>(pv);
    }

    static T * reallocate(T* old_ptr, const size_t n) {
        if (n == 0) {return nullptr;}
        if (n > static_cast<size_t>(-1) / sizeof(T)) {
            throw std::bad_array_new_length();
        }
        void * const pv = realloc(old_ptr, n * sizeof(T));
        if (!pv) { throw std::bad_alloc(); }
        return static_cast<T *>(pv);
    }

    void deallocate(T * const p, size_t) const noexcept {
        free(p);
#ifdef __linux__
        malloc_trim(0);
#endif
    }
};

#endif
