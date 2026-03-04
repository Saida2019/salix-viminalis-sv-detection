BEGIN { len=0 }
/^>/ {
    if(NR>1) print prev_chr "\t0\t" len
    split(substr($0,2), a, " ")
    prev_chr = a[1]
    len = 0
    next
}
{ len += length($0) }
END { print prev_chr "\t0\t" len }
