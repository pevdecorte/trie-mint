#ifndef TRIE_NODE_H_
#define TRIE_NODE_H_

#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

template <typename A>
class TrieNode {
  public:
    TrieNode(const std::string& path) :
        path_(path) {};
    TrieNode(const TrieNode&) = delete;
    TrieNode(TrieNode&&) = delete;
    void operator=(const TrieNode&) = delete;
    void operator=(TrieNode&&) = delete;

    // Returns a pointer to child c of this node; returns
    // nullptr if the child does not exist.
    TrieNode<A>* operator[](const char c) const {
        auto c_child = children_.find(c);
        if (c_child == children_.end()) {
            return nullptr;
        }
        return c_child->second.get();
    };

    // Adds the string held in array to the subtrie rooted
    // at this TrieNode. 
    // Nodes after the (terminal_after)th are marked as terminal nodes.
    void add_char_array(
            const char* array,
            const int terminal_after,
            const int terminal_until,
            const A value) {
        if (array[0] == '\0') return;

        // There won't be anymore terminal nodes, so we can
        // just bail out here.
        if (terminal_until == 0) return;

        if (children_.find(array[0]) == children_.end()) {
            children_[array[0]] = std::make_unique<TrieNode>(
                    path_ + array[0]);
            children_chars_.push_back(array[0]);
        }
        TrieNode* next_node = children_[array[0]].get();

        if (terminal_after <= 1) { 
            next_node->is_terminal_ = true;
            next_node->values_.push_back(value);
        }

        next_node->add_char_array(
                array + 1,
                terminal_after - 1,
                terminal_until - 1,
                value);
    };

    // Returns true iff this is a terminal node.
    bool is_terminal() const {
        return is_terminal_;
    };

    // Returns the path taken to this node.
    const std::string& path() const {
        return path_;
    };

    // Returns a pointer to the vector of this node's values.
    const std::vector<A>* values() const {
        return &values_;
    };

    // Returns a vector of all the char children.
    const std::vector<char>& children() const {
        return children_chars_;
    };

  private:
    bool is_terminal_ = false;
    std::unordered_map<char, std::unique_ptr<TrieNode>> children_;
    const std::string path_;  // Path to this node.
    std::vector<A> values_;
    std::vector<char> children_chars_;
};
#endif
