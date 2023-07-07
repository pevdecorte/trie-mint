// The binary accepts an input file consisting of patterns
// to match. It then reads a haystack string from the haystack
// file if one is provided on the command line, and otherwise
// from stdin.
// All matches are then printed to stdout.

/*
Each line of the patterns file should have the following format.

<pattern> <pattern name>



Each time a match is found, a line is printed to stdout having
the following format.

    <matched string> <start index of match in haystack>
        <start index of match in pattern> <name of pattern>

If the haystack file is provided on the command line, the
file name is echoed to stdout before the matches are printed.
*/


#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
#include <unordered_set>
#include <utility>

#include <stdio.h>

#include "gopt/gopt.h"
#include "matcher.h"
#include "trie.h"

struct Annotation {
    std::string name;
    int start;
};

// Reads the patterns file contents into the trie.
void ReadPatternsFileOrDie(
        Trie<Annotation>& trie,
        std::istream& patterns_file,
        const bool add_substrings = false,
        const int length_lower_bound = 0,
        const int length_upper_bound = 0);

// If the flag to ReadPatternsFileOrDie() is set indicating
// that all substrings within a given length range should
// also be added, this function does that.
void AddSubstringsIntoTrie(
        Trie<Annotation>& trie,
        const std::string& name,
        const std::string& s,
        int length_lower_bound,
        int length_upper_bound);

// Dies and reports that the patterns file was not properly
// formatted.
void DieOnInvalidFileFormat();

// Reads characters from the haystack and prints out
// match results.
// Returns true on success, false on failure.
bool ReadHaystack(
        std::istream* haystack,
        Matcher<Annotation>& matcher,
        const bool include_debug);

// Functions for displaying trie contents.
void PrintAllKeys(const Trie<Annotation>& patterns);
void PrintAllKeysAndValues(const Trie<Annotation>& patterns);

int main(int argc, char** argv) {
    struct option options[10];

    // Prints help information.
    options[0].long_name = "help";
    options[0].flags = GOPT_ARGUMENT_FORBIDDEN;

    // Sets the wildcard character for pattern strings.
    options[1].long_name = "pattern-wildcard";
    options[1].short_name = 'p';
    options[1].flags = GOPT_ARGUMENT_REQUIRED;

    // Sets the wildcard character for the haystack string.
    options[2].long_name = "haystack-wildcard";
    options[2].short_name = 'h';
    options[2].flags = GOPT_ARGUMENT_REQUIRED;

    // Prints out all the patterns and names from the patterns file
    // and then exits.
    options[3].long_name = "print-patterns-and-names";
    options[3].flags = GOPT_ARGUMENT_FORBIDDEN;

    // Prints out all the patterns from the patterns file and then exits.
    options[4].long_name = "print-patterns";
    options[4].flags = GOPT_ARGUMENT_FORBIDDEN;

    options[5].long_name = "range-lower";
    options[5].flags = GOPT_ARGUMENT_REQUIRED;

    options[6].long_name = "range-upper";
    options[6].flags = GOPT_ARGUMENT_REQUIRED;

    options[7].long_name = "suppress-header";
    options[7].flags = GOPT_ARGUMENT_FORBIDDEN;

    options[8].long_name = "no-debug";
    options[8].flags = GOPT_ARGUMENT_FORBIDDEN;

    options[9].flags = GOPT_LAST;

    argc = gopt(argv, options);
    gopt_errors(argv[0], options);

    std::ifstream patterns_file;
    std::unique_ptr<std::istream> haystack;

    // Include tRNA debugging labels in the output.
    const bool include_debug = !(options[8].count >= 1);


    if (options[0].count >= 1) {
        std::cout << "find_patterns <pattern file> [haystack file]"
            << " [--range-lower=... --range-upper=...]\n";
        return 0;
    }

    if (argc < 2) {
        std::cerr << "patterns file required\n";
        return 1;
    } else {
        patterns_file.open(argv[1]);
        if (!patterns_file.is_open()) {
            std::cerr << "Error opening patterns file " << argv[1]
                << std::endl;
            exit(1);
        }
    }

    if (argc >= 3) {
        haystack.reset(new std::ifstream(argv[2]));
        if (!haystack->good()) {
            std::cerr << "Error opening haystack file " << argv[2]
               << std::endl;
           exit(1);
        }
    } else {
        haystack.reset(&std::cin);
    } 


    if (options[5].count != options[6].count) {
        std::cerr << "Two range arguments must be provided,"
           << " or none at all.\n";
        exit(1);
    }

    // Parses range arguments from the command line.
    int length_lower_bound, length_upper_bound;
    bool length_bounds_set = false;
    if (options[5].count >= 1 && options[6].count >= 1) {
        std::istringstream iss1(options[5].argument);
        std::istringstream iss2(options[6].argument);
        iss1 >> length_lower_bound;
        iss2 >> length_upper_bound;
        if (iss1.fail() || iss2.fail()) {
            std::cerr << "Invalid range arguments.\n";
            exit(1);
        }
        length_bounds_set = true;
    }

    // Reads the patterns file and populates the Trie.
    Trie<Annotation> patterns;
    if (length_bounds_set) {
        ReadPatternsFileOrDie(patterns, patterns_file, true,
                length_lower_bound, length_upper_bound);
    } else {
        ReadPatternsFileOrDie(patterns, patterns_file);
    }
    patterns_file.close();

    // Prints all keys and values, then exits.
    if (options[3].count >= 1) {
        PrintAllKeysAndValues(patterns);
        exit(0);
    }

    // Prints all keys, then exits.
    if (options[4].count >= 1) {
        PrintAllKeys(patterns);
        exit(0);
    }

    // Constructs the Matcher object.
    Matcher<Annotation> matcher(patterns, include_debug);

    // Sets wildcard character if required.
    if (options[1].count >= 1) {
        std::cout << "Pattern wildcard " << options[1].argument
            << " provided. "
            << "Pattern wildcards are currently not supported.\n";
        exit(0);
    }
    if (options[2].count >= 1) {
        matcher.set_wildcard(options[2].argument[0]);
    }

    // Prints header if required.
    if (haystack.get() != &std::cin && options[7].count < 1) {
        std::cout << argv[2] << std::endl;
    }

    // Reads through haystack.
    if (!ReadHaystack(haystack.get(), matcher, include_debug)) {
        exit(1);
    }

    return 0;
}

bool ReadHaystack(
        std::istream* haystack,
        Matcher<Annotation>& matcher,
        const bool include_debug) {
    std::unordered_set<std::string> already_shown;

    char c;
    haystack->get(c);
    while (!haystack->fail()) {
        if (c == '\n') {
            haystack->get(c);
            continue;  // Skips over newlines.
        }
        matcher.advance(c);
        
        for (auto it = matcher.begin(); it != matcher.end(); ++it) {
            auto* match = it.match();
            std::string match_text;
            match_text += match->key + " " +
                std::to_string(match->start_index);

            if (include_debug) {
                match_text += " " + std::to_string(match->value.start) +
                    " " + match->value.name;
            }

            auto line_it = already_shown.find(match_text);
            if (line_it == already_shown.end()) {
                std::cout << match_text << std::endl;
                already_shown.insert(match_text);
            }
        }

        haystack->get(c);
    }

    if ((haystack->rdstate() & std::istream::eofbit) == 0) {
        std::cerr << "Unknown error encountered.\n";
        return false;
    }

    return true;
}

// Reads the patterns from patterns_file into the Trie provided.
// Each line in patterns_file is considered as a pattern.
// patterns_file should be a file open for input.
//
// This function may terminate the program with exit status 1
// if the patterns file cannot be parsed.
void ReadPatternsFileOrDie(
        Trie<Annotation>& trie,
        std::istream& patterns_file,
        const bool add_substrings,
        const int length_lower_bound,
        const int length_upper_bound) {
    std::string line;
    while (getline(patterns_file, line)) {
        std::istringstream iss(line);
        std::istream_iterator<std::string> eos;
        std::istream_iterator<std::string> field_it(iss);

        if (field_it == eos) { DieOnInvalidFileFormat(); }
        const std::string name = *field_it;

        ++field_it;
        if (field_it == eos) { DieOnInvalidFileFormat(); }
        const std::string pattern = *field_it;

        if (add_substrings) {
            AddSubstringsIntoTrie(
                    trie,
                    name,
                    pattern,
                    length_lower_bound,
                    length_upper_bound);
        } else {
            trie.add(pattern, Annotation{name});
        }
    }
}

void DieOnInvalidFileFormat() {
    std::cerr << "Invalid pattern file format\n";
    exit(1);
}

// Adds all substrings of lengths between length_lower_bound
// and length_upper_bound.
void AddSubstringsIntoTrie(
        Trie<Annotation>& trie,
        const std::string& name,
        const std::string& s,
        int length_lower_bound,
        int length_upper_bound) {
    const char* c_string = s.c_str();
    for (int i = 0; c_string[0] && i <= s.length() - length_lower_bound;
            ++i, ++c_string) {
        trie.add(c_string, length_lower_bound, length_upper_bound,
                Annotation{name, i});
    }
}

void PrintAllKeys(const Trie<Annotation>& patterns) {
    for (const std::string& s : patterns.AllKeys()) {
        std::cout << s << std::endl;
    }
}

void PrintAllKeysAndValues(const Trie<Annotation>& patterns) {
    for (const std::pair<std::string, Annotation>& s
            : patterns.AllKeysAndValues()) {
        std::cout << s.first << " " << s.second.name << std::endl;
    }
}

