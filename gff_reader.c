#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <limits.h>
#include "hash.h"
#include "gff_reader.h"

static int create_chromosome_transcript_hash(ChromosomeTranscript ***hash, unsigned long hash_size)
{
    *hash = malloc(sizeof(ChromosomeTranscript *) * hash_size);
    memset(*hash, 0, sizeof(ChromosomeTranscript *) * hash_size);
    return 0;
}

static int create_transcript_hash(Transcript ***hash, unsigned long hash_size)
{
    *hash = malloc(sizeof(Transcript *) * hash_size);
    memset(*hash, 0, sizeof(Transcript *) * hash_size);
    return 0;
}

static unsigned long add2transcript_hash(Transcript **hash, unsigned long hash_size, char *transcript, char *chromosome, char strand, unsigned long start, unsigned long end, char type)
{
    int flag = 0;
    unsigned long hash_value = ElfHash(transcript) % hash_size;
    Transcript *transcript_node = hash[hash_value];
    while (transcript_node)
    {
        if (!strcmp(transcript_node->transcript, transcript))
        {
            break;
        }
        transcript_node = transcript_node->next;
    }

    if (!transcript_node) // Add new transcript node.
    {
        flag = 1;
        transcript_node = malloc(sizeof(Transcript));
        transcript_node->transcript = malloc(sizeof(char) * (strlen(transcript) + 1));
        strcpy(transcript_node->transcript, transcript);
        transcript_node->chromosome = malloc(sizeof(char) * (strlen(chromosome) + 1));
        strcpy(transcript_node->chromosome, chromosome);
        transcript_node->strand = strand;
        transcript_node->start = 0;
        transcript_node->end = 0;
        transcript_node->exon_number = 0;
        transcript_node->cds_number = 0;
        transcript_node->element_number = 0;
        transcript_node->element = NULL;
        transcript_node->next = hash[hash_value];
        hash[hash_value] = transcript_node;
    }
    if (type == 'e')
    {
        transcript_node->exon_number++;
    }
    else
    {
        transcript_node->cds_number++;
    }
    Element *element_node = malloc(sizeof(Element));
    element_node->type = type;
    element_node->positions[0] = start;
    element_node->positions[1] = end;
    element_node->next = transcript_node->element;
    transcript_node->element = element_node;
    return flag;
}

static int add2chromosome_transcript_hash(ChromosomeTranscript **hash, unsigned long hash_size, char *chromosome, char *transcript, unsigned long plus)
{
    unsigned long hash_value = ElfHash(chromosome) % hash_size;
    ChromosomeTranscript *node = hash[hash_value];
    while (node)
    {
        if (!strcmp(node->chromosome, chromosome))
        {
            break;
        }
        node = node->next;
    }

    if (!node)
    {
        node = malloc(sizeof(ChromosomeTranscript));
        node->chromosome = malloc(sizeof(char) * strlen(chromosome) + 1);
        strcpy(node->chromosome, chromosome);
        node->transcript_number = 0;
        node->transcript_index = 0;
        node->transcript = NULL;
        node->next = hash[hash_value];
        hash[hash_value] = node;
    }
    node->transcript_number += plus;
    return 0;
}

ChromosomeTranscript *find_chromosome_transcript(ChromosomeTranscript **hash, unsigned long hash_size, char *key)
{
    unsigned long hash_value = ElfHash(key) % hash_size;
    ChromosomeTranscript *node = hash[hash_value];
    while (node)
    {
        if (!strcmp(node->chromosome, key))
        {
            break;
        }
        else
        {
            node = node->next;
        }
    }
    return node;
}

static int free_transcript_hash(Transcript **hash, unsigned long hash_size)
{
    for (unsigned long hash_index = 0; hash_index < hash_size; hash_index++)
    {
        for (Transcript *transcript_node = hash[hash_index]; transcript_node; transcript_node = transcript_node->next)
        {
            free(transcript_node->chromosome);
            free(transcript_node->transcript);
            Element *element_node = transcript_node->element;
            Element *element_node_ = NULL;
            while (element_node)
            {
                element_node_ = element_node->next;
                free(element_node);
                element_node = element_node_;
            }
        }
    }
    free(hash);
    return 0;
}

int free_chromosome_transcript_hash(ChromosomeTranscript **hash, unsigned long hash_size)
{
    for (unsigned long hash_index = 0; hash_index < hash_size; hash_index++)
    {
        for (ChromosomeTranscript *chromosome_transcript_node = hash[hash_index]; chromosome_transcript_node; chromosome_transcript_node = chromosome_transcript_node->next)
        {
            free(chromosome_transcript_node->chromosome);
            for (unsigned long transcript_index = 0; transcript_index < chromosome_transcript_node->transcript_number; transcript_index++)
            {
                free((chromosome_transcript_node->transcript + transcript_index)->transcript);
                free((chromosome_transcript_node->transcript + transcript_index)->element);
            }
        }
        if (hash[hash_index])
        {
            free(hash[hash_index]);
        }
    }
    free(hash);
    return 0;
}

static int compare_element(const void *x1, const void *x2)
{
    return ((Element *)x1)->positions[0] < ((Element *)x2)->positions[0] ? -1 : 1;
}

static int compare_transcript(const void *x1, const void *x2)
{
    return ((Transcript *)x1)->start < ((Transcript *)x2)->start ? -1 : 1;
}

static int parse_gene_structure(ChromosomeTranscript **chromosome_transcript_hash, Transcript **transcript_hash, unsigned long hash_size)
{
    unsigned long element_array_size = 0;
    for (unsigned long hash_index = 0; hash_index < hash_size; hash_index++)
    {
        /* Initialization. */
        for (ChromosomeTranscript *chromosome_transcript_node = chromosome_transcript_hash[hash_index]; chromosome_transcript_node; chromosome_transcript_node = chromosome_transcript_node->next)
        {
            chromosome_transcript_node->transcript = malloc(sizeof(Transcript) * chromosome_transcript_node->transcript_number);
        }
        /* #promoter + #exon + #intron + #cds + #5'utr + #3'utr */
        /* 1 + #exon + #exon - 1 + #exon + 2 */
        /* 3 * #exon + 2 */
        for (Transcript *transcript_node = transcript_hash[hash_index]; transcript_node; transcript_node = transcript_node->next)
        {
            if (element_array_size < transcript_node->exon_number)
            {
                element_array_size = transcript_node->exon_number;
            }
        }
    }
    element_array_size = 3 * element_array_size + 2;
    Element *element_array = malloc(sizeof(Element) * element_array_size);
    /*
    Structure of element_array:
        exon
        cds
        intron
        promoter if existed
        utr
    */
    for (unsigned long hash_index = 0; hash_index < hash_size; hash_index++)
    {
        for (Transcript *transcript_node = transcript_hash[hash_index]; transcript_node; transcript_node = transcript_node->next)
        {
            ChromosomeTranscript *chromosome_transcript_node = find_chromosome_transcript(chromosome_transcript_hash, hash_size, transcript_node->chromosome);
            Transcript *new_transcript_node = chromosome_transcript_node->transcript + (chromosome_transcript_node->transcript_index++);
            new_transcript_node->transcript = malloc(sizeof(char) * (strlen(transcript_node->transcript) + 1));
            strcpy(new_transcript_node->transcript, transcript_node->transcript);
            new_transcript_node->chromosome = NULL;
            new_transcript_node->strand = transcript_node->strand;
            new_transcript_node->exon_number = transcript_node->exon_number;
            new_transcript_node->cds_number = transcript_node->cds_number;
            new_transcript_node->element_number = 0;

            unsigned long cds_index = transcript_node->exon_number;

            for (Element *element = transcript_node->element; element; element = element->next)
            {
                if (element->type == 'e') // Add exon.
                {
                    memcpy(element_array + new_transcript_node->element_number, element, sizeof(Element));
                    (element_array + new_transcript_node->element_number)->next = NULL;
                    new_transcript_node->element_number++;
                }
                else // Add cds.
                {
                    memcpy(element_array + cds_index, element, sizeof(Element));
                    (element_array + cds_index)->next = NULL;
                    cds_index++;
                }
            }
            new_transcript_node->element_number += transcript_node->cds_number;

            if (transcript_node->exon_number)
            {
                qsort(element_array, transcript_node->exon_number, sizeof(Element), compare_element); // Sort exons.

                /* Add intron. */
                for (unsigned long element_index_ = 0; element_index_ < transcript_node->exon_number - 1; element_index_++)
                {
                    (element_array + new_transcript_node->element_number)->type = 'i';
                    (element_array + new_transcript_node->element_number)->next = NULL;
                    (element_array + new_transcript_node->element_number)->positions[0] = (element_array + element_index_)->positions[1] + 1;
                    (element_array + new_transcript_node->element_number)->positions[1] = (element_array + element_index_ + 1)->positions[0] - 1;
                    new_transcript_node->element_number++;
                }

                /* Add promoter. */
                if (transcript_node->strand == '-')
                {
                    (element_array + new_transcript_node->element_number)->positions[0] = (element_array + transcript_node->exon_number - 1)->positions[1] + 1;
                    (element_array + new_transcript_node->element_number)->positions[1] = (element_array + transcript_node->exon_number - 1)->positions[1] + PROMOTER;
                    (element_array + new_transcript_node->element_number)->type = 'p';
                    new_transcript_node->element_number++;
                }
                else if ((element_array + new_transcript_node->element_number)->positions[0] != 1)
                {
                    (element_array + new_transcript_node->element_number)->positions[0] = element_array->positions[0] > PROMOTER ? element_array->positions[0] - PROMOTER : 1;
                    (element_array + new_transcript_node->element_number)->positions[1] = element_array->positions[0] - 1;
                    (element_array + new_transcript_node->element_number)->type = 'p';
                    new_transcript_node->element_number++;
                }

                /* Add utr. */
                if (transcript_node->cds_number)
                {
                    unsigned long cds_min = ULONG_MAX;
                    unsigned long cds_max = 0;
                    for (unsigned long element_index_ = transcript_node->exon_number; element_index_ < transcript_node->exon_number + transcript_node->cds_number; element_index_++)
                    {
                        if (cds_min > (element_array + element_index_)->positions[0])
                        {
                            cds_min = (element_array + element_index_)->positions[0];
                        }
                        if (cds_max < (element_array + element_index_)->positions[1])
                        {
                            cds_max = (element_array + element_index_)->positions[1];
                        }
                    }

                    for (unsigned long element_index_ = 0; element_index_ < transcript_node->exon_number; element_index_++)
                    {
                        if ((element_array + element_index_)->positions[0] < cds_min)
                        {
                            (element_array + new_transcript_node->element_number)->type = (transcript_node->strand == '-' ? '3' : '5');
                            (element_array + new_transcript_node->element_number)->next = NULL;
                            (element_array + new_transcript_node->element_number)->positions[0] = (element_array + element_index_)->positions[0];
                            (element_array + new_transcript_node->element_number)->positions[1] = (element_array + element_index_)->positions[1] < cds_min ? (element_array + element_index_)->positions[1] : cds_min - 1;
                            new_transcript_node->element_number++;
                        }
                        if ((element_array + element_index_)->positions[1] > cds_max)
                        {
                            (element_array + new_transcript_node->element_number)->type = (transcript_node->strand == '-' ? '5' : '3');
                            (element_array + new_transcript_node->element_number)->next = NULL;
                            (element_array + new_transcript_node->element_number)->positions[0] = (element_array + element_index_)->positions[0] > cds_max ? (element_array + element_index_)->positions[1] : cds_max + 1;
                            (element_array + new_transcript_node->element_number)->positions[1] = (element_array + element_index_)->positions[1];
                            new_transcript_node->element_number++;
                        }
                    }
                }
            }
            qsort(element_array, new_transcript_node->element_number, sizeof(Element), compare_element);
            new_transcript_node->start = element_array->positions[0];
            new_transcript_node->end = (element_array + new_transcript_node->element_number - 1)->positions[1];
            new_transcript_node->element = malloc(sizeof(Element) * new_transcript_node->element_number);
            memcpy(new_transcript_node->element, element_array, sizeof(Element) * new_transcript_node->element_number);
            new_transcript_node->next = NULL;
        }
    }

    for (unsigned long hash_index = 0; hash_index < hash_size; hash_index++) // Sort.
    {
        for (ChromosomeTranscript *chromosome_transcript_node = chromosome_transcript_hash[hash_index]; chromosome_transcript_node; chromosome_transcript_node = chromosome_transcript_node->next)
        {
            qsort(chromosome_transcript_node->transcript, chromosome_transcript_node->transcript_number, sizeof(Transcript), compare_transcript);
        }
    }
    free(element_array);
    return 0;
}

int read_gff_file(char *file, ChromosomeTranscript ***chromosome_transcript_hash, unsigned long hash_size)
{
    char buffer[LINE];
    char ref_name[LINE]; // 0
    char type[LINE]; // 2
    unsigned long start; // 3
    unsigned long end; // 4
    char strand; // 6
    char attributes[LINE]; // 8

    create_chromosome_transcript_hash(chromosome_transcript_hash, hash_size);
    Transcript **transcript_hash = NULL;
    create_transcript_hash(&transcript_hash, hash_size);

    FILE *open_file = fopen(file, "r");
    while (fgets(buffer, LINE, open_file))
    {
        bool new_line = true;
        if ((buffer[0] != '#') && new_line)
        {
            const char *format = "%s\t%*[^\t]\t%s\t%lu\t%lu\t%*s\t%c\t%*c\t%[^\n]";
            sscanf(buffer, format, ref_name, type, &start, &end, &strand, attributes);
            if (!strcmp(type, "CDS") || !strcmp(type, "exon"))
            {
                char *attributes_pointer = attributes;
                char *sep = NULL;
                char transcript[LINE];
                while ((sep = strsep(&attributes_pointer, ";")))
                {
                    if(!strncmp(sep, "Parent=", 7))
                    {
                        strcpy(transcript, sep + 7);
                        break;
                    }
                }
                // printf("type: %s, transcript: %s, strand: %c, ref_name: %s, start: %lu, end: %lu, attributes: %s\n", type, transcript, strand, ref_name, start, end, attributes);
                unsigned long plus = add2transcript_hash(transcript_hash, hash_size, transcript, ref_name, strand, start, end, strcmp(type, "CDS") ? 'e' : 'c');
                add2chromosome_transcript_hash(*chromosome_transcript_hash, hash_size, ref_name, transcript, plus);
            }
        }
        new_line = buffer[strlen(buffer) - 1] == '\n' ? true : false;
    }
    fclose(open_file);
    parse_gene_structure(*chromosome_transcript_hash, transcript_hash, hash_size);
    free_transcript_hash(transcript_hash, hash_size);
    return 0;
}