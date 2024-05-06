trioPaternal=("NA19818" "HG01256" "HG00418" "NA19128" "NA12889" "NA12891")
trioMaternal=("NA19819" "HG01257" "HG00419" "NA19127" "NA12890" "NA12892")
trioChild=("NA19828" "HG01258" "HG00420" "NA19129" "NA12877" "NA12878")

map=$1
vcf=$2
out=$3

# on the filtered genotypes
for i in 0 1 2 3 4 5; do
    python scripts/check-mendelian-consistency.py -map $(map) -vcf $(vcf) -child ${trioChild[$i]} -father ${trioPaternal[$i]} -mother ${trioMaternal[$i]} -out $(out) > $(out)/${trioChild[$i]}.log;
done;