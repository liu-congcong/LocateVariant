CC = cc

LocateVariant: LocateVariant.c gff_reader.o gff_reader.h hash.o hash.h vcf_reader.o vcf_reader.h
	${CC} -Wall -g -o LocateVariant LocateVariant.c gff_reader.o hash.o vcf_reader.o

gff_reader.o: gff_reader.c gff_reader.h hash.o hash.h
	${CC} -Wall -c -o gff_reader.o gff_reader.c

hash.o: hash.c hash.h
	${CC} -Wall -c -o hash.o hash.c

vcf_reader.o: vcf_reader.c vcf_reader.h hash.o hash.h
	${CC} -Wall -c -o vcf_reader.o vcf_reader.c

clean:
	rm *.o