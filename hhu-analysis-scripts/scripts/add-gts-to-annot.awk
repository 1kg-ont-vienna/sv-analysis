#!/usr/bin/awk -f

BEGIN {
    FS="\t";
    OFS=FS;
    cnt1=0;
    cnt2=0;
}
{
    if (FNR == NR)
    {
        if ($1 ~ /^#/) 
        {
            next;
        }
        a[cnt1]=$3;
        cnt1 += 1;
    }
    else
    {
        if ($1 ~ /^#/)  
        {
            print $0;
            next;
        }
        if ($3==a[cnt2])
        {
            print $0;
            cnt2 += 1;
        }
    }
}