#
# 1. Preprocesses the chromosome .fa files and saves the resulting
# files to gen/chr.
#
# 2. Prepares the .patterns file with data from a .fa file.
#
# 3. Runs many-string-search/build/find_patterns for each of the
# chromosome .fa files in gen/chr. The matches found
# in <filename>.fa.mint are saved in gen/<filename>.matches.
#
# 4. Runs the postprocessing script which produces the mint lookup
# table and saves it in a .lookup file.
#
# All generated files are stored in gen/

NAMEROOT="rn6"
BASEDIR=${NAMEROOT}
DATAFILE="rn6-tRNAs-confidence-set.ss"
SWITCHES="-s"

set -x

# Makes sure the current build of find_patterns is up-to-date.
cd find_patterns
make
cd ..

# Preprocesses chromosome .fa files.
if [ ! -d gen/${BASEDIR} ]; then
    mkdir gen/${BASEDIR}
fi
if [ ! -d gen/${BASEDIR}/chr ]; then
    mkdir gen/${BASEDIR}/chr
fi
for f in $(ls data/${BASEDIR}/chr); do
    if [ -e gen/${BASEDIR}/chr/${f}.mint ]; then
        continue
    fi
    ./chr_preprocess.sh data/${BASEDIR}/chr/${f} \
        >gen/${BASEDIR}/chr/${f}.mint
done


# Prepares the .names file.
if [ ! -e gen/${BASEDIR}/${NAMEROOT}-tRNAs.names ]; then
    python3 naming.py ${SWITCHES} data/${BASEDIR}/${DATAFILE} |
        sort -k1,1 \
        >gen/${BASEDIR}/${NAMEROOT}-tRNAs.names
fi


# Prepares the .patterns file.
if [ ! -e gen/${BASEDIR}/${NAMEROOT}-tRNAs.patterns ]; then
    python3 patterns_and_intervals.py ${SWITCHES} \
        data/${BASEDIR}/${DATAFILE} \
        >gen/${BASEDIR}/${NAMEROOT}-tRNAs.patterns
fi

# Runs find_patterns on each preprocessed chromosome .mint file.
if [ ! -d gen/${BASEDIR}/matches ]; then
    mkdir gen/${BASEDIR}/matches
fi
for f in $(ls gen/${BASEDIR}/chr/*fa.mint); do
    filename=$(sed 's/.*\///' <<< ${f} | sed 's/\([^\.]*\).*/\1/')
    if [ -e gen/${BASEDIR}/matches/${filename}.matches ]; then
        continue
    fi  
    find_patterns/build/find_patterns \
        gen/${BASEDIR}/${NAMEROOT}-tRNAs.patterns \
        ${f} \
        --range-lower=16 \
        --range-upper=50 \
        >gen/${BASEDIR}/matches/${filename}.matches &
done

wait


# Runs reduced postprocessing to produce the .lookup file.
if [ ! -e gen/${BASEDIR}/${NAMEROOT}-tRNAs.lookup ]; then
    python3 reduced_postprocess.py ${SWITCHES} \
        data/${BASEDIR}/${DATAFILE} \
        gen/${BASEDIR}/${NAMEROOT}-tRNAs.names \
        gen/${BASEDIR}/matches/*.matches \
        >gen/${BASEDIR}/${NAMEROOT}-tRNAs.lookup
fi

