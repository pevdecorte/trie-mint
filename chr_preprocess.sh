# The preprocessing script is meant to be used with .fa
# chromosome files.
#
# It does the following:
#
# 
# 1. Removes the header line.
#
# 2. Uppercases all letters.
#
# 3. Removes newlines.
#
# 4. Prints the result to stdout.

tail -n +2 $1 | awk '{ print toupper($0) }' | tr -d '\n'
