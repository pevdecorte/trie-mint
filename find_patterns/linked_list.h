#ifndef LINKED_LIST_H_
#define LINKED_LIST_H_

#include <memory>

// Node class for linked list class.
// A pointer to the previous node is saved for fast deletion.

template <class T>
struct LinkedListNode {
    std::unique_ptr<T> value;
    LinkedListNode* previous;
    std::unique_ptr<LinkedListNode> next;
};

template <class T>
class LinkedListNodeIterator;

// Linked list class.
template <class T>
class LinkedList {
    friend class LinkedListNodeIterator<T>;
  public:
    // Creates a dummy node at the head.
    LinkedList() : head_(std::make_unique<LinkedListNode<T>>()),
        size_(0) {};

    LinkedList(const LinkedList&) = delete;
    void operator=(const LinkedList&) = delete;
    LinkedList(LinkedList&& rhs) :
        size_(rhs.size_), head_(std::move(rhs.head_)) {};
    void operator=(LinkedList&& rhs) {
        size_ = rhs.size_;
        head_ = std::move(rhs.head_);
    };


    // Adds a new node to the beginning of the linked list.
    void add(std::unique_ptr<T> new_value) {
        std::unique_ptr<LinkedListNode<T>> new_node =
            std::make_unique<LinkedListNode<T>>();
        new_node->value = std::move(new_value);
        new_node->previous = head_.get();

        new_node->next = std::move(head_->next);
        if (new_node->next != nullptr) {
            new_node->next->previous = new_node.get();
        }
        head_->next = std::move(new_node);
        ++size_;
    }

    void remove(LinkedListNode<T>* node) {
        if (node->next != nullptr) {
            node->next->previous = node->previous;
        }
        node->previous->next = std::move(node->next);
        --size_;
    }

    void remove(LinkedListNodeIterator<T>& node_iterator) {
        LinkedListNode<T>* previous = node_iterator.current_->previous;
        remove(node_iterator.current_);
        node_iterator.current_ = previous;
    }

    LinkedListNodeIterator<T> begin() const {
        return LinkedListNodeIterator<T>(this);
    }
    LinkedListNodeIterator<T> end() const {
        return LinkedListNodeIterator<T>();
    }

    int size() const {
        return size_;
    }

  private:
    std::unique_ptr<LinkedListNode<T>> head_;
    int size_;
};

// Iterator class for going through nodes in the linked
// list. Supports deletion of nodes while looping.
template <class T>
class LinkedListNodeIterator {
    friend class LinkedList<T>;
  public:
    LinkedListNodeIterator() :
        current_(nullptr) {}
    LinkedListNodeIterator(const LinkedList<T>* const list) :
        current_(list->head_->next.get()) {}

    bool operator!=(const LinkedListNodeIterator<T>& rhs) const {
        return !(current_ == rhs.current_);
    }

    LinkedListNodeIterator<T>& operator++() {
        if (current_ != nullptr) {
            current_ = current_->next.get();
        }
        return *this;
    }

    // Returns true if and only if we're not at the end of the
    // list.
    operator bool() const {
        return !(current_ == nullptr);
    }

    T* get_value() const {
        return current_->value.get();
    }

  private:
    LinkedListNode<T>* current_;
};

#endif
