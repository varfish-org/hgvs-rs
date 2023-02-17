BEGIN {
    OFS = FS;
    copy_table = "";

    # define list of known tables
    a = "gene|transcript|seq_anno|seq|associated_accessions|" \
        "exon-set|exon|exon_aln";
    split(a, names, "|")

    # read known ids
    for (i in names) {
        name = names[i];
        while ((getline line < ("download/" name ".tsv")) > 0) {
            known[name, line] = 1
            # print name, "|", line
        }
        close(name ".tsv")
    }
}

{
    if ($0 ~ /^COPY/) {
        gsub(/.*?\./, "", $2);
        table_name = $2;
        print $0;
    } else if ($0 ~ /^\\./) {
        table_name = "";
    }

    if (table_name == "" || known[table_name, $1] == 1) {
        print $0;
    }
}
