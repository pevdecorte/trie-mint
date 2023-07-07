#ifndef TRIE_H_
#define TRIE_H_

#include <memory>
#include <string>
#include <string.h>
#include <utility>
#include <vector>

#include "trie_node.h"

template <class A>
class Trie {
  public:
    Trie() : root_(std::make_unique<TrieNode<A>>("")) {} ;
    Trie(const Trie&) = delete;
    Trie(Trie&&) = delete;
    void operator=(const Trie&) = delete;
    void operator=(Trie&&) = delete;

    // Adds the string held in array to the trie, but counts
    // every node after the (terminal_after)th as a terminal node.
    void add(
            const char* array,
            const int terminal_after,
            const int terminal_until,
            const A value) {
        root_->add_char_array(array, terminal_after,
                terminal_until, value);
    };

    // Overloaded to handle C++ strings.
    void add(
            const std::string& str,
            const int terminal_after,
            const int terminal_until,
            const A value = A()) {
        add(str.c_str(), terminal_after, value);
    };

    void add(
            const std::string& str,
            const A value = A()) {
        add(str.c_str(), str.length(), str.length(), value);
    };

    // Returns the root node.
    const TrieNode<A>* root() const {
        return root_.get();
    };

    // Returns a vector containing all the strings that have
    // been added. Mainly for debugging.
    std::vector<std::string> AllKeys() const {
        std::vector<std::string> list;
        AllKeys(root_.get(), list);
        return list;
    };

    // Same as above but also includes values.
    std::vector<std::pair<std::string, A>> AllKeysAndValues() const {
        std::vector<std::pair<std::string, A>> list;
        AllKeysAndValues(root_.get(), list);
        return list;
    };

  private:
    const std::unique_ptr<TrieNode<A>> root_;

    void AllKeys(const TrieNode<A>* root,
            std::vector<std::string>& list) const {
        if (root->is_terminal()) {
            list.push_back(root->path());
        }

        for (const char& c : root->children()) {
            AllKeys((*root)[c], list);
        }
    };

    void AllKeysAndValues(
            const TrieNode<A>* root,
            std::vector<std::pair<
            std::string, A>>& list) const {
        if (root->is_terminal()) {
            for (const A& value : *root->values()) {
                list.emplace_back(make_pair(root->path(), value));
            }
        }

        for (const char& c : root->children()) {
            AllKeysAndValues((*root)[c], list);
        }
    };
};
#endif
