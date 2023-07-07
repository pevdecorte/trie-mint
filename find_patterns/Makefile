CC=g++

BUILDDIR = build
TESTS = linked_list_test trie_test slider_test \
        matcher_test wildcard_test
RUN_TESTS = $(patsubst %, %_run, $(TESTS))


$(BUILDDIR)/find_patterns: find_patterns.cc linked_list.h matcher.h \
    slider.h trie.h trie_node.h find_patterns.cc gopt/gopt.c \
    gopt/gopt-errors.c | $(BUILDDIR)
	$(CC) $< gopt/gopt.c gopt/gopt-errors.c -o $@

$(BUILDDIR):
	mkdir $(BUILDDIR)

$(BUILDDIR)/linked_list_test: linked_list_test.cc linked_list.h\
    slider.h trie.h trie_node.h | $(BUILDDIR)
	$(CC) $^ -o $@

$(BUILDDIR)/trie_test: trie_test.cc slider.h trie.h \
    trie_node.h | $(BUILDDIR)
	$(CC) $^ -o $@

$(BUILDDIR)/slider_test: slider_test.cc slider.h trie.h \
    trie_node.h | $(BUILDDIR)
	$(CC) $^ -o $@

$(BUILDDIR)/matcher_test: matcher_test.cc matcher.h linked_list.h \
    slider.h trie.h trie_node.h | $(BUILDDIR)
	$(CC) $^ -o $@

$(BUILDDIR)/wildcard_test: wildcard_test.cc matcher.h linked_list.h \
    slider.h trie.h trie_node.h | $(BUILDDIR)
	$(CC) $^ -o $@

$(RUN_TESTS): %_run: $(BUILDDIR)/%
	$^
	$^ | diff - test_output/$*.output

test: $(patsubst %, $(BUILDDIR)/%, $(TESTS)) 
	for t in $(TESTS); do \
            $(BUILDDIR)/$$t | diff - test_output/$$t.output; done

clean:
	rm -f -v $(BUILDDIR)/*.o $(BUILDDIR)/*test \
            $(BUILDDIR)/find_patterns
	rmdir build

