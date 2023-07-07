#include <iostream>
#include <string>

using namespace std;

#include "slider.h"
#include "trie.h"

class EmptyValue {};

int main() {
    Trie<EmptyValue> trie;
    string strings[] = {"the",
        "regal",
        "rainy",
        "rainiest",
        "rain",
        "in",
        "spain",
        "falls",
        "mainly",
        "in",
        "the",
        "plains"};

    for (auto s : strings) {
        cout << "Adding " << s << endl;
        trie.add(s);
    }
    cout << endl;

    for (auto s : trie.AllKeys()) {
        cout << s << endl;
    }
    return 0;
}
