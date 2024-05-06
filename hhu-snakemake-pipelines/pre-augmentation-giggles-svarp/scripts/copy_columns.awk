BEGIN {
    OFS=FS;
    count1=0;
    count2=0;
    print_info_header=0;
}
{
    # iterate over first vcf
    if (FNR==NR) {
        # skip header
        if (substr($0,1,1) == "#") {
            next;
        }
        #store column 1, 2, 3, and 8 of first vcf
        c1[count1]=$1;
        c2[count1]=$2;
        c3[count1]=$3;
        c8[count1]=$8;
        count1++;
        next;
    }
    # iterate over second vcf
    {
        if (substr($0,1,2) == "##") {
            print $0;
            next;
        }
        else if (substr($0,1,1) == "#") {
            print "##INFO=<ID=ID,Number=A,Type=String,Description=\"Variant IDs per ALT allele.\">";
            print $0;
            next;
        }
        # use assert to verify that the record is same.
        if (!(c1[count2]==$1)) {print "Column 1 match failed"; exit 1;}
        if (!(c2[count2]==$2)) {print "Column 2 match failed"; exit 1;}
        if (!(c3[count2]==$3)) {print "Column 3 match failed"; exit 1;}
        tmp=$8;
        # replace 8th column
        $8=c8[count2];
        # print whole line to stdout
        print $0;
        $8=tmp;
        count2++;
    }
}