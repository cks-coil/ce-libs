#include <vector>
#include <functional>
#include <iterator>


#ifndef __INCLUDED_LAZY_VECTOR_HPP__
#define __INCLUDED_LAZY_VECTOR_HPP__

template <typename _Tp> class LazyVectorIterator;
template <typename _Tp> class LazyVectorConstIterator;

template <typename _Tp>
class LazyVector: private std::vector< std::function<_Tp(void)> >{
public:
    typedef LazyVectorIterator<_Tp> iterator;
    typedef LazyVectorIterator<_Tp> const_iterator;
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
    _Tp at(int pos) const{ return (std::vector< std::function<_Tp(void)> >::at(pos))(); }
    _Tp operator[](int pos) const{ return (std::vector< std::function<_Tp(void)> >::operator[](pos))(); }
    _Tp front(void) const{ return (std::vector< std::function<_Tp(void)> >::front())(); }
    _Tp back(void) const{ return (std::vector< std::function<_Tp(void)> >::back())(); }
    iterator begin(void) { return iterator(this, 0); }
    iterator end(void) { return iterator(this, this->size()); }
    const_iterator begin(void) const { return const_iterator(this, 0); }
    const_iterator end(void) const { return const_iterator(this, this->size()); }
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
    _Tp operator*() const{
        return _lazyVector->at(_index);
    }
    bool operator==(const LazyVectorIterator<_Tp>& iterator) const{
        return (iterator._index == _index && iterator._lazyVector == _lazyVector);
    }
    bool operator!=(const LazyVectorIterator<_Tp>& iterator) const{
        return (iterator._index != _index || iterator._lazyVector != _lazyVector);
    }

private:
    LazyVectorIterator(LazyVector<_Tp> const *lazyVector, int index){
        _index = index;
        _lazyVector = lazyVector;
    }
    size_t _index;
    LazyVector<_Tp> const *_lazyVector;
};

#endif

