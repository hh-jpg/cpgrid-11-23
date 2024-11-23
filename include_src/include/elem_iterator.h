/*
 * @Author: jltu && jltu.siat.ac.cn, 1620808124@qq.com
 * @Date: 2024-07-16 09:42:58
 * @LastEditTime: 2024-08-28 10:16:17
 * @FilePath: /cpgrid/include/elem_iterator.h
 * @Description:
 *
 */
#ifndef ELEM_ITERATOR_H
#define ELEM_ITERATOR_H
#include <vector>
#include <iostream>
#include <iterator>
#include <functional>
#include "config.h"
// FilterIterator 类定义
template <typename Iterator, typename Predicate>
class FilterIterator {
public:
    using value_type = typename std::iterator_traits<Iterator>::value_type;
    using reference = typename std::iterator_traits<Iterator>::reference;
    using pointer = typename std::iterator_traits<Iterator>::pointer;
    using difference_type = typename std::iterator_traits<Iterator>::difference_type;
    using iterator_category = std::forward_iterator_tag;

    FilterIterator(Iterator begin, Iterator end, Predicate pred, int processor_id)
        : current(begin), end(end), predicate(pred), processor_id(processor_id) {
        advance_to_next_valid();
    }

    reference operator*() const {
        return *current;
    }

    pointer operator->() const {
        return &(*current);
    }

    FilterIterator& operator++() {
        ++current;
        advance_to_next_valid();
        return *this;
    }

    FilterIterator operator++(int) {
        FilterIterator temp = *this;
        ++(*this);
        return temp;
    }

    bool operator==(const FilterIterator &other) const {
        return current == other.current;
    }

    bool operator!=(const FilterIterator &other) const {
        return !(*this == other);
    }

private:
    void advance_to_next_valid() {
        while (current != end && !predicate(*current, processor_id)) {
            ++current;
        }
    }

    Iterator current;
    Iterator end;
    Predicate predicate;
    int processor_id;
};

// FilterRange 类定义
template<typename Iterator, typename Predicate>
class FilterRange {
public:
    FilterRange(Iterator begin, Iterator end, Predicate pred, int processor_id)
        : begin_(begin, end, pred, processor_id), end_(end, end, pred, processor_id) {}

    FilterIterator<Iterator, Predicate> begin() const {
        return begin_;
    }

    FilterIterator<Iterator, Predicate> end() const {
        return end_;
    }

private:
    FilterIterator<Iterator, Predicate> begin_;
    FilterIterator<Iterator, Predicate> end_;
};
#endif // ELEM_ITERATOR_H