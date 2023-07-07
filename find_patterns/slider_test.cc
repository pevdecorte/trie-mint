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

    cout << "Adding strings...\n";
    for (string s : strings) {
        trie.add(s);
        cout << s << endl;
    }
    cout << endl;

    cout << "children of root:\n";
    for (char c : trie.root()->children()) {
        cout << c << endl;
    }
    cout << endl;

    string host_string = "rainyspainymainlytheregalplainsyplain";

    Slider<EmptyValue> slider(trie.root(), 0);

    for (char c : host_string) {
        cout << c;
        if (!slider.feed(c)) {
            cout << " Slider is dead.\n";
            break;
        } else {
            cout << " matching...\n";
        }
        if (slider.has_match()) {
            cout << "Match found!\n";
            cout << "Matching string: \"" << slider.path()
                << "\"" << endl;
        }
    }

    return 0;
}
