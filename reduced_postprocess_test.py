# In theory, this should produce an identical file
# to gen/mm10/mm10-tRNAs.lookup.

python reduced_postprocess.py \
        data/mm10/mm10-tRNAs-confidence-set.ss \
        gen/mm10-reduced/mm10-reduced-tRNAs.names \
        gen/mm10-reduced/matches/*.matches \
        >gen/mm10-reduced/mm10-reduced-tRNAs.lookup-1

awk '{ print $1 " " $2 }' gen/mm10-reduced/mm10-reduced-tRNAs.lookup-1 \
        >gen/mm10-reduced/mm10-reduced-tRNAs.lookup

