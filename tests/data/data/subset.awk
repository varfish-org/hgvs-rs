BEGIN {
    OFS = FS;
    copy_table = "";

    # define list of known tables
    a = "gene|transcript|seq_anno|seq|associated_accessions|" \
        "exon_set|exon|exon_aln";
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
    if ($0 ~ /anonymous/) {
        next;  # skip for testing
    }

    # We don't need any of the materialized views.
    if ($0 ~ /REFRESH/ || $0 ~ /CREATE.*INDEX .*_mv/) {
        print "-- " $0;
        next;
    } else if ($0 ~ /MATERIALIZED VIEW/) {
        gsub(/MATERIALIZED VIEW/, "VIEW", $0);
        print $0;
        next;
    } else if ($0 ~ /WITH NO DATA/) {
        gsub(/WITH NO DATA/, "", $0);
        print $0;
        next;
    }

    if ($0 ~ /^COPY/) {
        table_name = $2;
        gsub(/.*?\./, "", table_name);
        if (table_name != "meta" && table_name != "origin") {
            print $0;  # do not print twice ;-)
        }
    } else if ($0 ~ /^\\./) {
        table_name = "";
    }

    if (table_name == "" || table_name == "meta" || \
        table_name == "origin" || known[table_name, $1] == 1) {
        print $0;
    }
}
