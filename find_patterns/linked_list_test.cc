#include <iostream>
#include <vector>

using namespace std;

#include "linked_list.h"
#include "slider.h"

class EmptyValue {};

int main() {
    LinkedList<Slider<EmptyValue>> list;

    for (int i = 0; i < 7; ++i) {
        list.add(std::make_unique<Slider<EmptyValue>>(nullptr, i));
    }
    for (auto item = list.begin(); item != list.end(); ++item) {
        cout << item.get_value()->start_index() << endl;
    }
    cout << endl;

    cout << "Delete everything.\n";
    for (auto it = list.begin(); it != list.end(); ++it) {
        list.remove(it);
    }

    cout << "Now add all the elements back.\n";
    for (int i = 0; i < 7; ++i) {
        list.add(std::make_unique<Slider<EmptyValue>>(nullptr, i));
    }
    cout << "Print all the elements.\n";
    for (auto item = list.begin(); item != list.end(); ++item) {
        cout << item.get_value()->start_index() << endl;
    }
    cout << endl;


    cout << "Delete 6 and 3.\n";
    for (auto item = list.begin(); item != list.end(); ++item) {
        cout << item.get_value()->start_index() << " ";
        if (item.get_value()->start_index() == 3) {
           list.remove(item);
        } else if (item.get_value()->start_index() == 6) {
            list.remove(item);
            continue;
        }
        cout << item.get_value()->start_index() << endl;
    }

    cout << endl;

    cout << "Print all the elements.\n";
    for (auto item = list.begin(); item != list.end(); ++item) {
        cout << item.get_value()->start_index() << endl;
    }
    
    cout << endl;

    cout << "Delete 0.\n";
    for (auto item = list.begin(); item != list.end(); ++item) {
        cout << item.get_value()->start_index() << " ";
        if (item.get_value()->start_index() == 0) {
           list.remove(item);
        }
        cout << item.get_value()->start_index() << endl;
    }

    cout << endl;

    cout << "Print all the elements.\n";
    for (auto item = list.begin(); item != list.end(); ++item) {
        cout << item.get_value()->start_index() << endl;
    }

    cout << endl;

    cout << "Delete 5.\n";
    for (auto item = list.begin(); item != list.end(); ++item) {
        if (item.get_value()->start_index() == 5) {
           list.remove(item);
        }
    }

    cout << endl;

    cout << "Print all the elements.\n";
    for (auto item = list.begin(); item != list.end(); ++item) {
        cout << item.get_value()->start_index() << endl;
    }

    return 0;
}
