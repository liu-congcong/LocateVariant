# LocateVariant

Locate variants in genes:

* Variant: variant id
* Chromosome: chromosome id
* Position: position (1-based)
* Transcript: transcript id | - (intergenic region)
* Type: exon | promoter | cds (**5'utr and 3'utr will be added soon**)
* Start: start position of the type (1-based)
* End: end position of the type (1-based)
* **Offset: offset will be added soon**

## Installation

```shell
git clone https://github.com/liu-congcong/LocateVariant.git
cd LocateVariant
gcc *.c -o LocateVariant
```

## Usage

```shell
LocateVariant -vcf VCF -gff GFF > OUTPUT
```
