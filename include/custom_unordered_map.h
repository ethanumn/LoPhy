#pragma once

#include <unordered_map>
#include <vector>
#include <string>

using std::unordered_map;
using std::string;
using std::vector;

/* Custom unordered_map class */
template<typename K, typename V>
class custom_unordered_map 
{
    private:
        unordered_map<K, V> map;

    public:
        V & operator[](const K & key) 
        {
            return map[key];
        }

        V & at(const K & key) 
        {
            return map.at(key);
        }

        void insert(const K & key, const V & value) 
        {
            map.insert({key, value});
        }

        template <typename InputIterator>
        void insert(InputIterator first, InputIterator last) 
        {
            map.insert(first, last);
        }

        bool has_key(const K & key) const 
        {
            return map.find(key) != map.end();
        }

        void erase(const K& key) 
        {
            map.erase(key);
        }

        auto begin() 
        {
            return map.begin();
        }

        auto end() 
        {
            return map.end();
        }

        auto begin() const 
        {
            return map.begin();
        }

        auto end() const 
        {
            return map.end();
        }

        size_t size() const 
        {
            return map.size();
        }

        bool empty() const 
        {
            return map.empty();
        }

        void clear() 
        {
            map.clear();
        }
};
