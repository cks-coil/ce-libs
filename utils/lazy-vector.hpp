#include <vector>
#include <functional>
#include <iterator>

#ifndef __INCLUDED_LAZY_VECTOR_HPP__
#define __INCLUDED_LAZY_VECTOR_HPP__

template <typename _Tp> class LazyVectorIterator;

template <typename _Tp>
class LazyVector: private std::vector< std::function<_Tp(void)> >{
public:
    typedef LazyVectorIterator<_Tp> iterator;
    using std::vector< std::function<_Tp(void)> >::assign;
    using std::vector< std::function<_Tp(void)> >::capacity;
    using std::vector< std::function<_Tp(void)> >::clear;
    using std::vector< std::function<_Tp(void)> >::get_allocator;
    using std::vector< std::function<_Tp(void)> >::insert;
    using std::vector< std::function<_Tp(void)> >::max_size;
    using std::vector< std::function<_Tp(void)> >::pop_back;
    using std::vector< std::function<_Tp(void)> >::push_back;
    using std::vector< std::function<_Tp(void)> >::reserve;
    using std::vector< std::function<_Tp(void)> >::resize;
    using std::vector< std::function<_Tp(void)> >::size;
    _Tp at(int pos){ return (std::vector< std::function<_Tp(void)> >::at(pos))(); }
    _Tp operator[](int pos){ return (std::vector< std::function<_Tp(void)> >::operator[](pos))(); }
    _Tp front(void){ return (std::vector< std::function<_Tp(void)> >::front())(); }
    _Tp back(void){ return (std::vector< std::function<_Tp(void)> >::back())(); }
    iterator begin(void){ return iterator(this, 0); }
    iterator end(void){ return iterator(this, this->size()); }
};

template <typename _Tp>
class LazyVectorIterator: public std::iterator<std::input_iterator_tag, _Tp>{
    friend LazyVector<_Tp>;
public:
    LazyVectorIterator(const LazyVectorIterator& iterator){
        _index = iterator._index;
        _lazyVector = iterator._lazyVector;
    }
    LazyVectorIterator<_Tp>& operator++(){
        _index++;
        return *this;
    }
    LazyVectorIterator<_Tp> operator++(int){
        LazyVectorIterator<_Tp>* result = *this;
        _index++;
        return result;
    }
    _Tp operator*(){
        return _lazyVector->at(_index);
    }
    bool operator==(const LazyVectorIterator<_Tp>& iterator){
        return (iterator._index == _index && iterator._lazyVector == _lazyVector);
    }
    bool operator!=(const LazyVectorIterator<_Tp>& iterator){
        return (iterator._index != _index || iterator._lazyVector != _lazyVector);
    }

private:
    LazyVectorIterator(LazyVector<_Tp> *lazyVector, int index){
        _index = index;
        _lazyVector = lazyVector;
    }
    size_t _index;
    LazyVector<_Tp> *_lazyVector;
};

#endif

