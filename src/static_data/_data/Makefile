ASSEMBLIES := \
	CHM1_1.0 \
	CHM1_1.1 \
	GRCh37 \
	GRCh37.p10 \
	GRCh37.p11 \
	GRCh37.p12 \
	GRCh37.p13 \
	GRCh37.p2 \
	GRCh37.p5 \
	GRCh37.p9 \
	GRCh38 \
	GRCh38.p10 \
	GRCh38.p11 \
	GRCh38.p12 \
	GRCh38.p1 \
	GRCh38.p2 \
	GRCh38.p3 \
	GRCh38.p4 \
	GRCh38.p5 \
	GRCh38.p6 \
	GRCh38.p7 \
	GRCh38.p8 \
	GRCh38.p9 \
	NCBI33 \
	NCBI34 \
	NCBI35 \
	NCBI36 \

JSONS := $(addsuffix .json.gz,${ASSEMBLIES})

.PHONY:
all: $(JSONS)

%.json.gz:
	wget -O $@ https://github.com/biocommons/bioutils/raw/main/src/bioutils/_data/assemblies/$@
