#ifndef SLIDER_H_
#define SLIDER_H_

// A "Slider" is effectively a pointer to a TrieNode. The Slider
// moves along the Trie as characters are fed to it, always
// pointing to a TrieNode matching all the characters fed to it
// so far.

#include <memory>
#include <string>
#include "trie.h"
#include "trie_node.h"

template <class A>
class Slider {
  public:
    Slider(const TrieNode<A>* start, const int start_index) :
        current_(start),
        start_index_(start_index) {};
    Slider(const Trie<A>& trie, const int start_index) :
        current_(trie.root()),
        start_index_(start_index) {};

    // Feeds one character to the slider, causing it to advance
    // inside the trie. Returns true upon success. If the currently
    // pointed to TrieNode has no c child, returns false.
    bool feed(const char c) {
        const TrieNode<A>* next_node = (*current_)[c];
        if (next_node == nullptr) {
            return false;
        }
        current_ = next_node;
        if (current_->is_terminal()) {
            match_found_ = true;
        } else {
            match_found_ = false;
        }
        ++match_length_;
        return true;
    };

    bool feed_wildcard(const char c) {
        ++wildcard_count_;
        ++match_length_;
        return feed(c);
    }

    // Returns true if a match has been found.
    bool has_match() const {
        return current_->is_terminal();
    };

    // Returns the path followed by the Slider. When has_match()
    // is true, this will be the value of the match found.
    const std::string& path() const {
        return current_->path();
    };

    const std::vector<A>* values() const {
        return current_->values();
    }

    int start_index() const { return start_index_; }

    const std::vector<char>& children() const {
        return current_->children();
    }

    int wildcard_count() const {
        return wildcard_count_;
    }

    int match_length() const {
        return match_length_;
    }

  private:
    // Currently pointed to node in a Trie.
    const TrieNode<A>* current_;
    // Position in haystack string where this slider started matching. 
    const int start_index_;
    // Are we currently sitting on the last character of a match?
    bool match_found_ = false;  
    // How many wildcards have we consumed so far?
    int wildcard_count_ = 0;
    // How many characters have been matched?
    int match_length_ = 0;
};

#endif
