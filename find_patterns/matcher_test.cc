#include <iostream>
#include <string>

using namespace std;

#include "matcher.h"
#include "slider.h"
#include "trie.h"

struct StringValue {
    std::string name = "";
};

int main() {
    Trie<StringValue> trie;
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
    int count = 1;
    for (string s : strings) {
        trie.add(s, StringValue{std::to_string(count)});
        ++count;
        cout << s << endl;
    }
    cout << endl;

    string host_string = "rainyspainymainlytheregalplainsyplain";
    cout << "Host string: " << host_string << endl;

    Matcher<StringValue> matcher(trie);

    for (char c : host_string) {
        matcher.advance(c);
        for (MatchIterator<StringValue> it = matcher.begin();
                it != matcher.end(); ++it) {
            const Match<StringValue>* match = it.match();
            cout << "Match found: " << match->key << " at index "
                << match->start_index << "; name: " << match->value.name
                << "; length: " << match->length
                << endl;
        }
    }

    return 0;
}
