#ifndef __GFF_READER_H__

#define __GFF_READER_H__

#define LINE 10240

#define FILE_NAME 10240

#define PROMOTER 2000

typedef struct Element
{
    char type;
    unsigned long positions[2];
    struct Element *next;
} Element;

typedef struct Transcript
{
    char *transcript;
    char *chromosome;
    char strand;
    unsigned long start;
    unsigned long end;
    unsigned long exon_number;
    unsigned long cds_number;
    unsigned long element_number;
    Element *element;
    struct Transcript *next;
} Transcript;

typedef struct ChromosomeTranscript
{
    char *chromosome;
    unsigned long transcript_number;
    unsigned long transcript_index;
    Transcript *transcript;
    struct ChromosomeTranscript *next;
} ChromosomeTranscript;

int read_gff_file(char *, ChromosomeTranscript ***, unsigned long);

ChromosomeTranscript *find_chromosome_transcript(ChromosomeTranscript **, unsigned long, char *);

int free_chromosome_transcript_hash(ChromosomeTranscript **, unsigned long);

#endif