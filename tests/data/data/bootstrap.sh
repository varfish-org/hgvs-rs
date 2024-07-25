#!/usr/bin/bash

# Setup Logging -------------------------------------------------------------

log()
{
    >&2 echo $@
}

debug()
{
    [[ "${VERBOSE-0}" -ne 0 ]] && >&2 echo $@
}

set -euo pipefail

if [[ "${VERBOSE-0}" -ne 0 ]]; then
    set -x
fi

psql-uta()
{
    echo "set schema '$VERSION'; $1" \
    | PGPASSWORD=anonymous psql --csv -h uta.biocommons.org -U anonymous -d uta \
    | tail -n +3 \
    | sort
}

pg-list()
{
    echo $* \
    | tr ' ' '\n' \
    | sed -e "s/^/'/g" -e "s/$/'/g" \
    | tr '\n' ',' \
    | sed -e 's/,$//g' \
    | sed -e "s/^/(/g" -e "s/$/)/g"
}

# Initialization ------------------------------------------------------------

if [[ "$#" -ne 2 ]]; then
    log "USAGE: bootstrap.sh DL_URL VERSION"
    log ""
    log "E.g.: bootstrap.sh http://dl.biocommons.org/uta uta_20210129"
    log ""
    log "Set VERBOSE=1 to increase verbosity."
    exit 1
fi

# path to the directory where the script resides.
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# Download URL.
DL_URL=$1

# Database name/version.
VERSION=$2

# Destination directory.
DST=$SCRIPT_DIR

# The HGNC symbols of the genes to fetch.
set +e
read -r -d '' GENES <<EOF
AADACL3
ADGRL3
ADRA2B
ADRB2
AGBL5
AKR1C3
ALG9
ALMS1
AOAH
ASB18
ATM
BAHCC1
BOLA3
BRCA1
BRCA2
C18orf32
CACNA1S
CCNQ
CEP120
COL1A1
COL3A1
DEFB133
DMKN
DNAH11
DNMBP-AS1
EGLN1
EGR2
FAM102A
FAM129B
FAM58A
FBN1
FBXO22
FOLR3
GTF3C2
H2AW
HELQ
HEXB
HIST3H2A
HMGA1
IFNA16
IRAK3
JRK
KCNIP4
KRAS
L1CAM
LCE2B
LCE3C
LINC00336
LPHN3
MC1R
MECP2
MED21
MLH1
MSH6
MYH7
NEFL
NR1H2
OMA1
OPA1
OR9A4
ORAI1
PLOD1
POMC
PTEN
PTH2
RABGAP1L
RBM38
RET
RGL3
RMRP
RPGR
RPL17
RPL38
RPL41
RYR1
RYR2
SCN5A
SDHC
SDR5A2
SELENON
SELL
SERPINC1
SIL1
SLC22A5
SLC34A3
SPG11
SRA1
SRD5A2
SSPO
SSTR3
TAZ
TERC
TPM3
TSC2
TSIX
UFD1
VPS26A
ZCCHC3
EOF
set -e

# Transform gene list for postgres query.
PG_GENES=$(pg-list $GENES)

# Force consistent sorting.
export LC_ALL=C

# Generate Data for Tests ---------------------------------------------------

cd $DST

# download database dumps

mkdir -p download
cd download
for f in $VERSION.pgd.gz{,.sha1}; do
    test -e $f || wget $DL_URL/$f
done
cd ..

# extract identifiers

psql-uta "select hgnc from gene g where g.hgnc in $PG_GENES;" >download/gene.tsv
psql-uta "select ac from transcript where hgnc in $PG_GENES;" >download/transcript.tsv
psql-uta "select seq_anno_id from seq_anno where ac in (select ac from transcript where hgnc in $PG_GENES);" >download/seq_anno.tsv
psql-uta "select seq_id from seq where seq_id in (select seq_id from seq_anno where ac in (select ac from transcript where hgnc in $PG_GENES));" >download/seq.tsv
psql-uta "select associated_accession_id from associated_accessions where tx_ac in (select ac from transcript where hgnc in $PG_GENES);" >download/associated_accessions.tsv
psql-uta "select exon_set_id from exon_set where tx_ac in (select ac from transcript where hgnc in $PG_GENES);" >download/exon_set.tsv
psql-uta "select exon_id from exon where exon_set_id in (select exon_set_id from exon_set where tx_ac in (select ac from transcript where hgnc in $PG_GENES));" >download/exon.tsv
psql-uta "select exon_aln_id from exon_aln where tx_exon_id in (select exon_id from exon where exon_set_id in (select exon_set_id from exon_set where tx_ac in (select ac from transcript where hgnc in $PG_GENES)));" >download/exon_aln.tsv

# build sql subset

pigz -d -c download/$VERSION.pgd.gz \
| awk -F ' ' -f subset.awk \
| pigz -c \
> $VERSION-subset.pgd.gz
