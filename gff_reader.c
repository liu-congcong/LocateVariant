#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
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
    unsigned long flag = 1;
    unsigned long hash_value = ElfHash(transcript) % hash_size;
    Transcript *transcript_node = hash[hash_value];
    while (transcript_node)
    {
        if (!strcmp(transcript_node->transcript, transcript))
        {
            flag = 0; // Exist.
            break;
        }
        else if (transcript_node->next)
        {
            transcript_node = transcript_node->next;
        }
        else
        {
            break;
        };
    };
    if (flag)
    {
        Transcript *new_transcript_node = malloc(sizeof(Transcript));
        new_transcript_node->transcript = malloc(sizeof(char) * (strlen(transcript) + 1));
        strcpy(new_transcript_node->transcript, transcript);
        new_transcript_node->chromosome = malloc(sizeof(char) * (strlen(chromosome) + 1));
        strcpy(new_transcript_node->chromosome, chromosome);
        new_transcript_node->strand = strand;
        new_transcript_node->start = 0;
        new_transcript_node->end = 0;
        new_transcript_node->exon_number = 0;
        new_transcript_node->cds_number = 0;
        new_transcript_node->element = NULL;
        new_transcript_node->next = NULL;
        if (transcript_node)
        {
            transcript_node->next = new_transcript_node;
        }
        else
        {
            hash[hash_value] = new_transcript_node;
        };
        transcript_node = new_transcript_node;
    };
    if (type == 'e')
    {
        transcript_node->exon_number++;
    }
    else
    {
        transcript_node->cds_number++;
    };
    Element *new_element_node = malloc(sizeof(Element));
    new_element_node->type = type;
    new_element_node->positions[0] = start;
    new_element_node->positions[1] = end;
    new_element_node->next = NULL;
    Element *element_node = transcript_node->element;
    if (element_node)
    {
        while (element_node->next)
        {
            element_node = element_node->next;
        };
        element_node->next = new_element_node;
    }
    else
    {
        transcript_node->element = new_element_node;
    };
    return flag;
}

static int add2chromosome_transcript_hash(ChromosomeTranscript **hash, unsigned long hash_size, char *chromosome, char *transcript, unsigned long plus)
{
    int flag = false;
    unsigned long hash_value = ElfHash(chromosome) % hash_size;
    ChromosomeTranscript *node = hash[hash_value];
    while (node)
    {
        if (!strcmp(node->chromosome, chromosome))
        {
            flag = true;
            break;
        }
        else if (node->next)
        {
            node = node->next;
        }
        else
        {
            break;
        };
    };
    if (!flag)
    {
        ChromosomeTranscript *new_node = malloc(sizeof(ChromosomeTranscript));
        new_node->chromosome = malloc(sizeof(char) * strlen(chromosome) + 1);
        strcpy(new_node->chromosome, chromosome);
        new_node->transcript_number = 0;
        new_node->transcript_index = 0;
        new_node->transcript = NULL;
        new_node->next = NULL;
        if (node)
        {
            node->next = new_node;
        }
        else
        {
            hash[hash_value] = new_node;
        };
        node = new_node;
    };
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
        };
    };
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
            };
        };
    };
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
            };
        };
        if (hash[hash_index])
        {
            free(hash[hash_index]);
        };
    };
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
    for (unsigned long hash_index = 0; hash_index < hash_size; hash_index++) // Initialization.
    {
        for (ChromosomeTranscript *chromosome_transcript_node = chromosome_transcript_hash[hash_index]; chromosome_transcript_node; chromosome_transcript_node = chromosome_transcript_node->next)
        {
            chromosome_transcript_node->transcript = malloc(sizeof(Transcript) * chromosome_transcript_node->transcript_number);
        };
    };

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

            unsigned long element_array_size = transcript_node->exon_number * 2 + transcript_node->cds_number;
            Element *element_array = malloc(sizeof(Element) * element_array_size); // #exon + #cds + #intron + #promoter = 2 * #exon + #cds
            unsigned long left_index = 0;
            unsigned long right_index = element_array_size - 1;
            if (transcript_node->exon_number)
            {
                for (Element *element = transcript_node->element; element; element = element->next) // Step1: add all exons to generate introns.
                {
                    if (element->type == 'e')
                    {
                        memcpy(element_array + left_index, element, sizeof(Element));
                        (element_array + left_index)->next = NULL;
                        left_index++;
                    }
                    else
                    {
                        memcpy(element_array + right_index, element, sizeof(Element));
                        (element_array + right_index)->next = NULL;
                        right_index--;
                    };
                };
                qsort(element_array, left_index, sizeof(Element), compare_element); // All exons have been sorted.

                if (transcript_node->strand == '-') // Step2: add promoter.
                {
                    (element_array + right_index)->positions[0] = (element_array + left_index - 1)->positions[1] + 1;
                    (element_array + right_index)->positions[1] = (element_array + left_index - 1)->positions[1] + 2000;
                }
                else
                {
                    (element_array + right_index)->positions[0] = element_array->positions[0] >= 2000 ? element_array->positions[0] - 2000 : 0;
                    (element_array + right_index)->positions[1] = element_array->positions[0] - 1;
                };
                (element_array + right_index)->type = 'p';
                right_index--;

                for (unsigned long left_index_ = 0; left_index_ < left_index - 1; left_index_++) // Add introns
                {
                    (element_array + right_index)->type = 'i';
                    (element_array + right_index)->next = NULL;
                    (element_array + right_index)->positions[0] = (element_array + left_index_)->positions[1] + 1;
                    (element_array + right_index)->positions[1] = (element_array + left_index_ + 1)->positions[0] - 1;
                    right_index--;
                };
                
            };
            qsort(element_array, element_array_size, sizeof(Element), compare_element);
            new_transcript_node->start = element_array->positions[0];
            new_transcript_node->end = (element_array + element_array_size - 1)->positions[1];
            new_transcript_node->element = element_array;
            new_transcript_node->next = NULL;
        };
    };

    for (unsigned long hash_index = 0; hash_index < hash_size; hash_index++) // Sort.
    {
        for (ChromosomeTranscript *chromosome_transcript_node = chromosome_transcript_hash[hash_index]; chromosome_transcript_node; chromosome_transcript_node = chromosome_transcript_node->next)
        {
            qsort(chromosome_transcript_node->transcript, chromosome_transcript_node->transcript_number, sizeof(Transcript), compare_transcript);
        };
    };
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
                    };
                };
                // printf("type: %s, transcript: %s, strand: %c, ref_name: %s, start: %lu, end: %lu, attributes: %s\n", type, transcript, strand, ref_name, start, end, attributes);
                unsigned long plus = add2transcript_hash(transcript_hash, hash_size, transcript, ref_name, strand, start, end, strcmp(type, "CDS") ? 'e' : 'c');
                add2chromosome_transcript_hash(*chromosome_transcript_hash, hash_size, ref_name, transcript, plus);
            };
        };
        new_line = buffer[strlen(buffer) - 1] == '\n' ? true : false;
    };
    fclose(open_file);
    parse_gene_structure(*chromosome_transcript_hash, transcript_hash, hash_size);
    free_transcript_hash(transcript_hash, hash_size);
    return 0;
}