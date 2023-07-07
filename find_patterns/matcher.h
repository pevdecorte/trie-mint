#ifndef MATCHER_H_
#define MATCHER_H_

#include <iterator>
#include <memory>
#include <string>

// DEBUG
#include <iostream>

#include "linked_list.h"
#include "slider.h"
#include "trie.h"
#include "trie_node.h"

template <class A>
struct Match {
    const std::string& key;
    const A& value;
    const int start_index;
    const int length;
};

// Wrapper class for LinkedListNodeIterator<Slider> which skips
// over Sliders that aren't currently sitting on a match.
template <class A>
class MatchIterator {
  public:
    MatchIterator() {};
    MatchIterator(const LinkedList<Slider<A>>* list, const bool include_debug) :
        include_debug_(include_debug) {
        trie_node_values_ = nullptr;
        // Searches for the first Slider with a match.
        for (current_ = list->begin(); current_ != list->end();
                ++current_) {
            if (current_.get_value()->has_match()) {
                trie_node_values_ = current_.get_value()->values();
                vector_it_ = trie_node_values_->begin();
                update_current_match();
                return;
            }
        }
    }

    const Slider<A>* slider() { return current_.get_value(); }

    bool operator!=(const MatchIterator<A>& rhs) const {
        return current_ != rhs.current_;
    }

    // If there are still values stored at the current TrieNode,
    // moves the pointer to the next one.
    // Otherwise, advances the current iterator to the next slider
    // sitting on a match.
    MatchIterator<A>& operator++() {
        if (include_debug_) {
            if (trie_node_values_ != nullptr) {
                ++vector_it_;
            }
            if (vector_it_ != trie_node_values_->end()) {
                update_current_match();
                return *this;
            }
            trie_node_values_ = nullptr;
        }
        ++current_;
        for (; current_ && !current_.get_value()->has_match();
                ++current_) {}
        if (current_) {
            trie_node_values_ = current_.get_value()->values();
            vector_it_ = trie_node_values_->begin();
            update_current_match();
        }
        return *this;
    }
    operator bool() const {
        return (bool)current_;
    } 

    const Match<A>* match() const {
        if (!current_) {
            return nullptr;
        }
        return current_match_.get();
    }

  private:
    // current_ is always either set to a node containing a
    // match, or it evaluates to false.
    LinkedListNodeIterator<Slider<A>> current_;
    const std::vector<A>* trie_node_values_;
    typename std::vector<A>::const_iterator vector_it_;
    std::unique_ptr<Match<A>> current_match_;
    const bool include_debug_ = true;

    void update_current_match() {
        current_match_ = std::unique_ptr<Match<A>>(
            new Match<A>{
                current_.get_value()->path(),
                *vector_it_,
                current_.get_value()->start_index(),
                current_.get_value()->match_length()
            }
        );
    };
};

// This class holds the trie, and keeps track of the current
// position in the string being searched.
template <class A>
class Matcher {
    friend class MatchIterator<A>;

  public:
    Matcher(const Trie<A>& fragments, const bool include_debug = true) :
        fragments_(fragments), start_index_(0), include_debug_(include_debug) {};
    Matcher() = delete;
    Matcher(const Matcher<A>&) = delete;
    Matcher(Matcher<A>&&) = delete;
    void operator=(const Matcher<A>&) = delete;
    void operator=(Matcher<A>&&) = delete;

    MatchIterator<A> begin() const {
        return MatchIterator<A>(&sliders_, include_debug_);
    }

    MatchIterator<A> end() const {
        return MatchIterator<A>();
    } 

    void advance(const char c) {
        sliders_.add(std::make_unique<Slider<A>>(
                    fragments_, start_index_));
        ++start_index_;
        if (has_wildcard_ && c == wildcard_) {
            LinkedList<Slider<A>> new_list;

            // If there have been too many wildcards in a row,
            // we simply kill all existing sliders.
            ++wildcards_in_a_row_;
            if (wildcard_limit_ != -1 &&
                    wildcards_in_a_row_ > wildcard_limit_) {
                sliders_ = std::move(new_list);
                return;
            }

            for (auto slider_it = sliders_.begin();
                    slider_it != sliders_.end();
                    ++slider_it) {
                for (const char c : slider_it.get_value()->children()) {
                    std::unique_ptr<Slider<A>> new_slider(
                            new Slider<A>(*slider_it.get_value()));
                    new_slider->feed_wildcard(c);
                    new_list.add(std::move(new_slider));
                }
            }
            sliders_ = std::move(new_list);
        } else {
            wildcards_in_a_row_ = 0;
            for (auto slider_it = sliders_.begin();
                    slider_it != sliders_.end();
                    ++slider_it) {
                if (!slider_it.get_value()->feed(c)) {
                    sliders_.remove(slider_it);
                }
            }
        }
    }

    int size() const {
        return sliders_.size();
    }

    void set_wildcard(const char c) {
        wildcard_ = c;
        has_wildcard_ = true;
    }

    void set_wildcard_limit(const int limit) {
        wildcard_limit_ = limit;
    }

  private:
    const Trie<A>& fragments_;
    LinkedList<Slider<A>> sliders_;
    char wildcard_;
    bool has_wildcard_;
    int wildcard_limit_ = -1;
    int wildcards_in_a_row_ = 0;
    int start_index_;
    const bool include_debug_;
};

#endif 
