#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

// CONSTANTS
#define MAX_KMER_LEN 100
#define MAX_LINE_LEN 1024
#define MAX_ADDED_LEN 12
#define SHIFTED_GAP_POSITIONS 0    // count also gap positions that are shifted from center
#define MAX_GAP_LENGTH 10

/* GLOBAL VARIABLES */
uint64_t mask_ULL[42][42];   /* PRIMARY mask_ULL FOR EACH SEPARATE NUCLEOTIDE */
uint64_t lowmask_ULL[42];    /* LOW MASK FOR EACH KMER */
uint64_t highmask_ULL[42];   /* HIGH MASK FOR EACH KMER */
short int Nlength = 32;              /* LENGTH OF MASK */

/* STRUCTURE DEFINITIONS */
typedef struct {
    char *sequence;
    int kmer_length;
    int background_count;
    int signal_count;
} LocalMaxEntry;

typedef struct {
    LocalMaxEntry *entries;
    int num_entries;
} LocalMaxTable;

typedef struct {
    int original_kmer_number;
    int original_kmer_value;
    int original_kmer_length;
    int original_gap_position;
    int original_gap_length;
    int added_length;
    int max_length;
    int triplegap_kmer_value;
    int gap_pos[3];
    int gap_len[3];
    int debug;
} SingleAndTripleGapKmer;

typedef struct {
   __uint128_t kmer_value;
   short int kmer_length;
   short int gap_position;
   short int gap_length;
} MergedKmer;

/* SUBROUTINE THAT PRINTS KMER (USES NO MASK) */
void Kmerprint (__uint128_t print_sequence_value, short int kmer_length, short int gap_position, short int gap_length)
{
    for (kmer_length--; kmer_length >= 0; kmer_length--, gap_position--)
    {   if (gap_position == 0) for(int i = 0; i < gap_length; i++) printf("n");
        printf("%c","ACGT"[(print_sequence_value >> (kmer_length * 2)) & 3]);
    }
}

void Triplegap_kmerprint(SingleAndTripleGapKmer* triplegapkmer) {
    // Calculate gap lengths based on original kmer
    (*triplegapkmer).gap_len[0] = (*triplegapkmer).original_gap_position;
    (*triplegapkmer).gap_len[1] = (*triplegapkmer).original_gap_length + (*triplegapkmer).gap_pos[0] - (*triplegapkmer).gap_pos[2];
    (*triplegapkmer).gap_len[2] = (*triplegapkmer).original_kmer_length + (*triplegapkmer).original_gap_length -
                                  (*triplegapkmer).gap_len[0] - (*triplegapkmer).gap_len[1] +
                                  (*triplegapkmer).gap_pos[0] - (*triplegapkmer).gap_pos[2];

    // Print padding for alignment
    for (int i = 0; i <= (*triplegapkmer).max_length - (*triplegapkmer).gap_pos[0]; i++) {
        printf(" ");
    }

    // Iterate through the positions in the triplegap kmer
    for (int i = 0; i <= (*triplegapkmer).added_length; i++) {
        // Print the first gap section (dots)
        if (i == (*triplegapkmer).gap_pos[0]) {
            for (int j = 0; j < (*triplegapkmer).gap_len[0]; j++) {
                printf(".");
            }
        }

        // Print the second gap section (n's for original gap content)
        if (i == (*triplegapkmer).gap_pos[1]) {
            for (int j = 0; j < (*triplegapkmer).gap_len[1]; j++) {
                printf("n");
            }
        }

        // Print the third gap section (dots)
        if (i == (*triplegapkmer).gap_pos[2]) {
            for (int j = 0; j < (*triplegapkmer).gap_len[2]; j++) {
                printf(".");
            }
        }

        // Print the bases from the triplegap kmer
        if (i < (*triplegapkmer).added_length) {
            // Extract the correct base and print
            int base = ((*triplegapkmer).triplegap_kmer_value >> (((*triplegapkmer).added_length - i - 1) * 2)) & 3;
            printf("%c", "ACGT"[base]);
        }
    }

    // Flush the output to ensure the print is displayed
    // fflush(stdout);
}

/* SUBROUTINE THAT GENERATES 128 bit BITMASKS FOR NUCLEOTIDES, KMER STRINGS AND DELETIONS */
void GenerateMask ()
{
    short int counter;
    short int start_position;
    short int position;
    short int current_kmer_length;
    /* PRIMARY mask_ULL FOR EACH SEPARATE NUCLEOTIDE */
    for (mask_ULL[1][0] = 3ULL, counter = 0; counter < Nlength; counter++) mask_ULL[1][counter+1] = mask_ULL[1][counter] << 2;
    /* GENERATES mask_ULLS FOR EXTRACTION OF EACH KMER STRING */
    for(current_kmer_length = 2; current_kmer_length <= Nlength; current_kmer_length++)
    {
        for (start_position = 0; start_position < Nlength-current_kmer_length; start_position++)
        {
            for (mask_ULL[current_kmer_length][start_position]=mask_ULL[1][start_position], position = start_position+1; position < current_kmer_length+start_position; position++) {mask_ULL[current_kmer_length][start_position] += mask_ULL[1][position];}
        }
    }
    /* GENERATES HIGH AND LOW mask_ULLS FOR DELETIONS */
    for (lowmask_ULL[0] = mask_ULL[1][0], position = 1; position < Nlength-2; position++) lowmask_ULL[position] = lowmask_ULL[position-1]+mask_ULL[1][position];
    for (highmask_ULL[Nlength-2] = mask_ULL[1][Nlength-2], position = Nlength-3; position > 0; position--) highmask_ULL[position] = highmask_ULL[position+1]+mask_ULL[1][position];
}

// Binary search in sorted array
int find_index(long int *indices, int n, long int target) {
   int l = 0, r = n - 1;
   while (l <= r) {
       int m = (l + r) / 2;
       if (indices[m] == target) return 1;
       if (indices[m] < target) l = m + 1;
       else r = m - 1;
   }
   return 0;
}

// Merges original kmer and triplegap kmer (in SingleAndTripleGapKmer structure) to one merged kmer
void Merge_kmers(SingleAndTripleGapKmer* triplegap_kmer, MergedKmer* merged, short int print_kmer) {
    // Initialize merged kmer value
    merged->kmer_value = 0;
    merged->kmer_length = triplegap_kmer->original_kmer_length + triplegap_kmer->added_length;

    int original_pos = 0; // Position tracker for the original kmer
    int current_pos = 0;  // Position tracker for the merged kmer

    for (int i = 0; i <= triplegap_kmer->added_length; i++) {
        // Handle first gap section (dots)
        if (i == triplegap_kmer->gap_pos[0]) {
            for (int j = 0; j < triplegap_kmer->gap_len[0]; j++) {
                merged->kmer_value = (merged->kmer_value << 2) |
                    ((triplegap_kmer->original_kmer_value >> ((triplegap_kmer->original_kmer_length - original_pos - 1) * 2)) & 3);
                original_pos++;
                current_pos++;
            }
        }

        // Handle middle gap section (n's)
        if (i == triplegap_kmer->gap_pos[1]) {
            merged->gap_position = current_pos;
            merged->gap_length = triplegap_kmer->gap_len[1];
            //if (merged->gap_length == 0) merged->gap_position = 0;  // Sets gap position to 0 if gap length is 0 to make representation compatible with rest of the code
            current_pos += triplegap_kmer->gap_len[1];
        }

        // Handle second gap section (dots)
        if (i == triplegap_kmer->gap_pos[2]) {
            for (int j = 0; j < triplegap_kmer->gap_len[2]; j++) {
                merged->kmer_value = (merged->kmer_value << 2) |
                    ((triplegap_kmer->original_kmer_value >> ((triplegap_kmer->original_kmer_length - original_pos - 1) * 2)) & 3);
                original_pos++;
                current_pos++;
            }
        }

        // Insert bases from the triplegap kmer
        if (i < triplegap_kmer->added_length) {
            merged->kmer_value = (merged->kmer_value << 2) |
                ((triplegap_kmer->triplegap_kmer_value >> ((triplegap_kmer->added_length - i - 1) * 2)) & 3);
            current_pos++;
        }
    }

    // Handle remaining original kmer bases after gaps
    while (original_pos < triplegap_kmer->original_kmer_length) {
        merged->kmer_value = (merged->kmer_value << 2) |
            ((triplegap_kmer->original_kmer_value >> ((triplegap_kmer->original_kmer_length - original_pos - 1) * 2)) & 3);
        original_pos++;
        current_pos++;
    }
    //if (merged->gap_length == 0) printf (" %i ", merged->gap_position);
    
    // Print the merged kmer for verification
    if (print_kmer == 1) Kmerprint(merged->kmer_value, merged->kmer_length, merged->gap_position, merged->gap_length);
}

// Helper function to extract a base from kmer_value
int Get_Base(__uint128_t kmer, int pos, int length) {
    // Each base is represented by 2 bits
    // Position 0 is the leftmost base
    return (kmer >> ((length - pos - 1) * 2)) & 3;
}


// Convert a position in the read containing the gap to a position in the bit representation of the gapped kmer
int Map_read_pos_to_gapped_kmer(int pos, int gap_pos, int gap_len) {
    if (pos <= gap_pos || gap_len == 0 || gap_pos == 0) return pos;
    else return pos - gap_len;
}

// Compare bases at given positions in two kmers
int Compare_bases(long long kmer1, int pos1, int len1,
                 long long kmer2, int pos2, int len2) {
    return Get_Base(kmer1, pos1, len1) == Get_Base(kmer2, pos2, len2);
}

// Helper function to check if position is in gap
int Is_in_gap(int pos, int gap_pos, int gap_len) {
    return (pos > gap_pos && pos <= gap_pos + gap_len);
}

// Compare bases in gapped kmers using actual positions from mapper
int Compare_bases_in_gapped_kmers(long long kmer1, int pos1, int gap_pos1, int gap_len1,
                                 long long kmer2, int pos2, int gap_pos2, int gap_len2,
                                 int len1, int len2) {
    
    // Jump over gap in original kmer by adjusting pos 1 and 2 if pos 1 is over gap pos 1
    if (pos1 > gap_pos1)
    {
        pos1 += gap_len1;
        pos2 += gap_len1;
    }
    
    /*/ Check if positions are in gaps
    int in_gap1 = Is_in_gap(pos1, gap_pos1, gap_len1);
    int in_gap2 = Is_in_gap(pos2, gap_pos2, gap_len2);
    
    if (in_gap1) return 1;  // N matches anything
    if (in_gap2) return in_gap1;  // X only matches N */
    
    // Remap to get actual positions in bit representation
    int actual_pos1 = Map_read_pos_to_gapped_kmer(pos1, gap_pos1, gap_len1);
    int actual_pos2 = Map_read_pos_to_gapped_kmer(pos2, gap_pos2, gap_len2);
    
    return Compare_bases(kmer1, actual_pos1, len1, kmer2, actual_pos2, len2);
}

// kmer aligner
int old_align_kmers(long long original_kmer, int original_length, int original_gap_pos, int original_gap_len,
              long long merged_kmer, int merged_length, int merged_gap_pos, int merged_gap_len) {
   
   // Try each possible starting position in merged kmer
   for (int start_pos = 0; start_pos <= merged_length - original_length; start_pos++) {
       int match = 1;
       
       // For each position in original kmer space
       for (int i = 0; i < original_length; i++) {
           // Compare using gapped positions, adjusting merged position by start_pos
           if (!Compare_bases_in_gapped_kmers(
                   original_kmer, i, original_gap_pos, original_gap_len,
                   merged_kmer, i + start_pos,
                   merged_gap_pos, merged_gap_len,
                   original_length, merged_length)) {
               match = 0;
               break;
           }
       }
       
       if (match) return start_pos;
   }
   
   return -1;
}

// kmer aligner
int align_kmers(long long original_kmer, int original_length, int original_gap_pos, int original_gap_len,
                long long merged_kmer, int merged_length, int merged_gap_pos, int merged_gap_len) {
    // Try each possible starting position in merged kmer
    int start_pos;
    int match = 1;
    for (start_pos = 0 ; start_pos <= merged_length - original_length; start_pos++) {
        //printf("%i:", start_pos);
        match = 1;
        int offset = start_pos;
        // For each position in original kmer space
        for (int i = 0; i < original_length; i++) {
            //printf("%i-", i);
            if (i == original_gap_pos) {offset += original_gap_len - merged_gap_len; /*printf("*");*/}
            if (!Compare_bases(original_kmer, i, original_length,
                               merged_kmer, offset + i, merged_length)) {match--; break;}
        }
        if (match) return start_pos;
    }
    return -1;
}

// helper function to map read position to gapped kmer in MergedKmer structure
int Map_read_pos_to_merged_kmer(MergedKmer* m, int read_pos)
{
    if (read_pos <= m->gap_position || m->gap_length == 0 || m->gap_position == 0) return read_pos;
    else return read_pos-m->gap_length;
}

// Function to extract triplegap sequence value from MergedKmer
long long Get_triplegap_sequence_value(SingleAndTripleGapKmer* o, MergedKmer* m) {
   long long triplegap_value = 0;
   int read_pos = 0;
   
   // Before first gap
   for (int i = 0; i < o->gap_pos[0]; i++) {
       int base = Get_Base(m->kmer_value, Map_read_pos_to_merged_kmer(m, read_pos++), m->kmer_length);
       triplegap_value = (triplegap_value << 2) | base;
   }
   read_pos += o->gap_len[0];
   
   // Between first and second gap
   for (int i = o->gap_pos[0]; i < o->gap_pos[1]; i++) {
       int base = Get_Base(m->kmer_value, Map_read_pos_to_merged_kmer(m, read_pos++), m->kmer_length);
       triplegap_value = (triplegap_value << 2) | base;
   }
   read_pos += o->gap_len[1];
   
   // Between second and third gap
   for (int i = o->gap_pos[1]; i < o->gap_pos[2]; i++) {
       int base = Get_Base(m->kmer_value, Map_read_pos_to_merged_kmer(m, read_pos++), m->kmer_length);
       triplegap_value = (triplegap_value << 2) | base;
   }
   read_pos += o->gap_len[2];
   
   // After third gap
   for (int i = o->gap_pos[2]; i < o->added_length; i++) {
       int base = Get_Base(m->kmer_value, Map_read_pos_to_merged_kmer(m, read_pos++), m->kmer_length);
       triplegap_value = (triplegap_value << 2) | base;
   }
   
   return triplegap_value;
}

// one-sided kmer alignment routine where gap of original kmer is N that can match any base or X, and gap in merged kmer is X, which can match only N
int Align_to_identify_triplegapkmer(SingleAndTripleGapKmer* o, MergedKmer* m) {
   // Find alignment position
    
    if(o->debug == 1) {
        printf (" *** aligning ");
        Kmerprint(o->original_kmer_value, o->original_kmer_length, o->original_gap_position, o->original_gap_length);
        printf(" to ");
        Kmerprint(m->kmer_value, m->kmer_length, m->gap_position, m->gap_length);
    }
   int pos = align_kmers(o->original_kmer_value, o->original_kmer_length,
                        o->original_gap_position, o->original_gap_length,
                        m->kmer_value, m->kmer_length,
                        m->gap_position, m->gap_length);
   
    if(o->debug == 1) printf (" align pos: %i ", pos);
    if (pos < 0)
    {
        if(o->debug == 1) printf ("\t**FAILED TO ALIGN** ");
        return 0;
    }
        
    short int original_right_half_length = o->original_kmer_length - o->original_gap_position;
    short int fill = o->original_gap_length - m->gap_length;
    short int left_fill = m->gap_position - o->original_gap_position - pos;
    if (m->gap_length == 0) left_fill = 0;  // if the gap is completely filled, m->gap position is zero which would mess up the calculation
    short int right_fill = fill - left_fill;
    short int total_length_difference = m->kmer_length + m->gap_length - o->original_kmer_length - o->original_gap_length;
    
    short int gap_position_start = o->original_gap_position + pos;  // position in merged kmer coordinate where original gap starts
    
    o->added_length = m->kmer_length - o->original_kmer_length;
    
    o->gap_pos[0] = pos;
    o->gap_pos[1] = o->gap_pos[0] + left_fill;
    o->gap_pos[2] = o->gap_pos[0] + fill;
    

    //o->gap_pos[2] = o->gap_pos[1] + right_fill;
    //o->gap_pos[2] = o->added_length - (total_length_difference - pos);
    //o->gap_pos[2] = m->gap_length - o->original_gap_length - pos;  // math from above
    //o->gap_pos[2] = o->added_length - (o->gap_pos[0] - fill + left_fill);
    

    o->gap_len[0] = o->original_gap_position;
    //o->gap_len[1] = m->gap_length - left_fill;
    o->gap_len[1] = o->original_gap_length - fill;
    o->gap_len[2] = original_right_half_length;
   
    //printf("\nCalculated triplegap positions: %d,%d,%d lengths: %d,%d,%d",o->gap_pos[0], o->gap_pos[1], o->gap_pos[2],o->gap_len[0], o->gap_len[1], o->gap_len[2]); fflush(stdout);
    
    //printf (" %i ", o->triplegap_kmer_value);
    //Kmerprint(m->kmer_value, m->kmer_length, m->gap_position, m->gap_length);
    o->triplegap_kmer_value = Get_triplegap_sequence_value(o, m);
    //printf ("-%i ", o->triplegap_kmer_value);
    
    if(o->debug > 1)
    {
        // Test by merging back
        MergedKmer merged_result;
        Merge_kmers(o, &merged_result, 0);
        
        printf("\nMerge test - Original: ");
        Triplegap_kmerprint(o);
        printf("\n");
        Kmerprint(m->kmer_value, m->kmer_length, m->gap_position, m->gap_length);
        printf(" Result: ");
        Kmerprint(merged_result.kmer_value, merged_result.kmer_length, merged_result.gap_position, merged_result.gap_length);
        fflush(stdout);
        
        //if (merged_result.kmer_value != m->kmer_value || merged_result.kmer_length != m->kmer_length || merged_result.gap_position != m->gap_position || merged_result.gap_length != m->gap_length) {return 0;}
    }
    return 1;
}

// Subroutine to find correct triplegap kmer when merged kmer and original kmer are known. Used for index mapping from 5d kmer index using file number, kmer length, gap position, gap length and kmer value to triplegap style 7d array using file number, original kmer number, triplegap kmer length, the three gap positions and triplegap kmer value.
short int Identify_triplegapkmer(SingleAndTripleGapKmer* o, MergedKmer* m) {
   
    // printf("\nKmer ");
    // Kmerprint (m->kmer_value, m->kmer_length, m->gap_position, m->gap_length);
    // printf("  Index %i,%i,%i converted to ", m->kmer_length, m->gap_position, m->gap_length);
    // Define triplegap_length, gap_fill, bases_added_to_ends
    int triplegap_length = m->kmer_length - o->original_kmer_length;
    o->added_length = triplegap_length;
    int gap_fill = o->original_gap_length - m->gap_length;
    int total_append = o->added_length - gap_fill;

    
    MergedKmer merged_result;
    // printf("Input triplegap Gap Positions: %d, %d, %d\n", o->gap_pos[0], o->gap_pos[1], o->gap_pos[2]);
    // printf("Input triplegap Gap Lengths: %d, %d, %d\n", o->gap_len[0], o->gap_len[1], o->gap_len[2]);
    
    // printf("\n Merged kmer: ");
    // Kmerprint(m->kmer_value, m->kmer_length, m->gap_position, m->gap_length);
    
    // printf("\tTriplegap Length: %d\t", triplegap_length);
    // printf("Gap Fill: %d\t", gap_fill);
    // printf("Total Bases Added to Ends: %d\n", total_append);

    // Iterate through all possible triplegap shifts
    for(int left_append = total_append; left_append >= 0; left_append--) {
        
        int invalid_pattern = 0;
        // Define gap positions and lengths
        o->gap_pos[0] = left_append;
        o->gap_len[0] = o->original_gap_position;
        
        o->gap_pos[1] = m->gap_position - o->original_gap_position;
        int left_side_gap_fill_length = o->gap_pos[1] - o->gap_pos[0];
        // disallows if gap positions not ordered
        if (left_side_gap_fill_length < 0) invalid_pattern = 10;
        // disallows if too many bases added to gap
        if (left_side_gap_fill_length > gap_fill) invalid_pattern++;
        o->gap_len[1] = o->original_gap_length - gap_fill + left_side_gap_fill_length;
        
        o->gap_pos[2] = o->added_length - left_append - 1;
        // disallows if wrong number of bases added to gap
        if (o->gap_pos[2] - o->gap_pos[0] != gap_fill) invalid_pattern += 2;
        o->gap_len[2] = o->original_kmer_length - o->original_gap_position;

        // printf("\n--- Testing triplegap with %d bases on the left ---\n", left_append);
        // printf("Triplegap Gap Positions: %d, %d, %d\n", o->gap_pos[0], o->gap_pos[1], o->gap_pos[2]);
        // printf("Triplegap Gap Lengths: %d, %d, %d\n", o->gap_len[0], o->gap_len[1], o->gap_len[2]);
        
        // if (invalid_pattern) printf ("\nInvalid pattern %i", invalid_pattern);
        
        // Extract triplegap_kmer_value from merged_kmer_value based on gap_pos and gap_len
        // printf("\n");
        long long triplegap_kmer_value = 0;
        int current_pos = 0;

        // Extract bases before first gap
        for (int i = 0; i < o->gap_pos[0]; i++) {
            if (current_pos >= m->kmer_length) {
                // printf("Invalid position %d. Skipping.\n", current_pos);
                triplegap_kmer_value = -1;
                break;
            }
            int base = Get_Base(m->kmer_value, current_pos, m->kmer_length);
            triplegap_kmer_value = (triplegap_kmer_value << 2) | base;
            // printf("%c", "ACGT"[base]);
            current_pos++;
        }
        if (triplegap_kmer_value == -1) continue;

        // Skip first gap
        current_pos += o->gap_len[0];
        // printf("-");

        // Extract bases between first and second gap
        for (int i = o->gap_pos[0]; i < o->gap_pos[1]; i++) {
            if (current_pos >= m->kmer_length) {
                // printf("Invalid position %d. Skipping.\n", current_pos);
                triplegap_kmer_value = -1;
                break;
            }
            int base = Get_Base(m->kmer_value, current_pos, m->kmer_length);
            triplegap_kmer_value = (triplegap_kmer_value << 2) | base;
            // printf("%c", "ACGT"[base]);
            current_pos++;
        }
        if (triplegap_kmer_value == -1) continue;

        // Skip second gap
        current_pos += o->gap_len[1];
        // printf("-");

        // Extract bases between second and third gap
        for (int i = o->gap_pos[1]; i < o->gap_pos[2]; i++) {
            if (current_pos >= m->kmer_length) {
                // printf("Invalid position %d. Skipping.\n", current_pos);
                triplegap_kmer_value = -1;
                break;
            }
            int base = Get_Base(m->kmer_value, current_pos, m->kmer_length);
            triplegap_kmer_value = (triplegap_kmer_value << 2) | base;
            // printf("%c", "ACGT"[base]);
            current_pos++;
        }
        if (triplegap_kmer_value == -1) continue;

        // Skip third gap
        current_pos += o->gap_len[2];
        // printf("-");

        // Extract bases after third gap
        for (int i = o->gap_pos[2]; i < o->added_length; i++) {
           if (current_pos >= m->kmer_length) {
               triplegap_kmer_value = -1;
               break;
           }
           int base = Get_Base(m->kmer_value, current_pos, m->kmer_length);
           triplegap_kmer_value = (triplegap_kmer_value << 2) | base;
           current_pos++;
        }
        if (triplegap_kmer_value == -1) continue;
        
        // Define a triplegap_kmer SingleAndTripleGapKmer struct
        SingleAndTripleGapKmer triplegap_kmer;
        triplegap_kmer.original_kmer_value = o->original_kmer_value;
        triplegap_kmer.original_kmer_length = o->original_kmer_length;
        triplegap_kmer.original_gap_position = o->original_gap_position;
        triplegap_kmer.original_gap_length = o->original_gap_length;
        triplegap_kmer.added_length = triplegap_length;
        triplegap_kmer.max_length = o->max_length;
        triplegap_kmer.triplegap_kmer_value = triplegap_kmer_value;
        triplegap_kmer.gap_pos[0] = o->gap_pos[0];
        triplegap_kmer.gap_pos[1] = o->gap_pos[1];
        triplegap_kmer.gap_pos[2] = o->gap_pos[2];
        triplegap_kmer.gap_len[0] = o->gap_len[0];
        triplegap_kmer.gap_len[1] = o->gap_len[1];
        triplegap_kmer.gap_len[2] = o->gap_len[2];

        //printf("\n");
        //Triplegap_kmerprint(&triplegap_kmer);
        //printf("  for triplegap Kmer Value: %d\n", triplegap_kmer.triplegap_kmer_value);

        // Call Merge_kmers
        Merge_kmers(&triplegap_kmer, &merged_result, 0);

        //Kmerprint(merged_result.kmer_value, merged_result.kmer_length, merged_result.gap_position, merged_result.gap_length);
        
        // printf("\tMerged Kmer Value: %llu\t", (unsigned long long)merged_result.kmer_value);
        // printf("Length: %d\t", merged_result.kmer_length);
        // printf("Gap Position: %d\t", merged_result.gap_position);
        // printf("Gap Length: %d\n", merged_result.gap_length);

        // Compare to input merged kmer
        if(merged_result.kmer_value == m->kmer_value &&
           merged_result.kmer_length == m->kmer_length &&
           merged_result.gap_position == m->gap_position &&
           merged_result.gap_length == m->gap_length)
        {
            //printf("\nMatch found when %d bases on the left!\n", left_append);
            // Print the triplegap kmer
            //Triplegap_kmerprint(&triplegap_kmer);
            //printf("**match** at %i,%i,%i,%i", o->added_length, o->gap_pos[0], o->gap_pos[1], o->gap_pos[2]);
            //printf("-- because of match %lu:%lu,%i:%i,%i:%i,%i:%i", (long int) merged_result.kmer_value, (long int) m->kmer_value, merged_result.kmer_length, m->kmer_length, merged_result.gap_position, m->gap_position, merged_result.gap_length, m->gap_length);
            return 1;
        }
    }

    //printf("\nNo matching triplegap kmer found.\n");
    //printf("NO MATCH at len %i:%i,%i,%i", o->added_length, o->gap_pos[0], o->gap_pos[1], o->gap_pos[2]);
    //printf("-- because of mismatch %lu:%lu,%i:%i,%i:%i,%i:%i", (long int) merged_result.kmer_value, (long int) m->kmer_value, merged_result.kmer_length, m->kmer_length, merged_result.gap_position, m->gap_position, merged_result.gap_length, m->gap_length);
    //printf(" Correct: ");
    //Kmerprint (m->kmer_value, m->kmer_length, m->gap_position, m->gap_length);
    //printf(" vs ");
    //Kmerprint (merged_result.kmer_value, merged_result.kmer_length, merged_result.gap_position, merged_result.gap_length);
    return 0;
}


// Index packer to generate one int describing the four indexing parameters
long int pack_index(int l, int len, int p1, int p2, int p3) {
    // Each field is shifted enough to avoid overlap, assuming they stay within expected ranges
    return (long int)l   * 1000000000000L
         + (long int)len * 1000000000L
         + (long int)p1  * 100000L
         + (long int)p2  * 1000L
         + (long int)p3;
}

// check if memory allocated. This does not work because alignment can generate sets of gap positions that are not allocated
int Kmer_out_of_range(long int *******counts, SingleAndTripleGapKmer *gap_iter, int file_number)
{
    // Check top-level pointer for file
    if (!counts || !counts[file_number]) return 1;

    // Check second dimension: original_kmer_number
    if (!counts[file_number][gap_iter->original_kmer_number]) return 1;

    // Check third dimension: added_length
    if (!counts[file_number][gap_iter->original_kmer_number][gap_iter->added_length]) return 1;

    // Check gap positions
    if (!counts[file_number][gap_iter->original_kmer_number][gap_iter->added_length][gap_iter->gap_pos[0]]) return 1;
    if (!counts[file_number][gap_iter->original_kmer_number][gap_iter->added_length][gap_iter->gap_pos[0]][gap_iter->gap_pos[1]]) return 1;
    if (!counts[file_number][gap_iter->original_kmer_number][gap_iter->added_length][gap_iter->gap_pos[0]][gap_iter->gap_pos[1]][gap_iter->gap_pos[2]]) return 1;

    // If we get here, it's allocated
    return 0;
}

int Kmer_out_of_range_binary(long int *******counts, SingleAndTripleGapKmer *gap_iter, int file_number) {
   long int *index_array = counts[0][0][0][0][0][0];
   if (!index_array) return 1;
   
   long int idx = pack_index(gap_iter->original_kmer_number, gap_iter->added_length, gap_iter->gap_pos[0],
                           gap_iter->gap_pos[1], gap_iter->gap_pos[2]);
   
   int left = 2;
   int right = (int)index_array[1];
   printf("Search idx %ld in array size %ld, ", idx, index_array[1]); fflush(stdout);
   
   while (left <= right) {
       int mid = (left + right) / 2;
       long int mid_idx = index_array[mid];
       
       if (mid_idx == idx) {
           printf("FOUND at %d\n", mid); fflush(stdout);
           return 0;
       }
       if (mid_idx < idx) left = mid + 1;
       else right = mid - 1;
   }
   printf("not found\n"); fflush(stdout);
   return 1;
}

// simple linear search for counted triplegap kmers (to prevent access of unallocated memory)
int Kmer_out_of_range_linear(long int *******counts, SingleAndTripleGapKmer *gap_iter, int file_number) {
   long int *index_array = counts[0][0][0][0][0][0];
   if (!index_array) return 1;
   
   long int idx = pack_index(gap_iter->original_kmer_number, gap_iter->added_length, gap_iter->gap_pos[0],
                           gap_iter->gap_pos[1], gap_iter->gap_pos[2]);
   
    
   //printf("Linear search idx %ld in allocation size %ld, array size %ld, compared to index ", idx, index_array[0], index_array[1]); fflush(stdout);
   
// number of objects is: counts[0][0][0][0][0][0][1] - 1
   for (int i = 2; i <= index_array[1]; i++) {
       //printf ("%ld ", index_array[i]);
       if (index_array[i] == idx) {
         //  printf("FOUND at %d\n", i); fflush(stdout);
           return 0;
       }
   }
   //printf("not found\n"); fflush(stdout);
   return 1;
}


// Comparison function for qsort (descending order)
static int compare_counts(const void *a, const void *b) {
    return (*(long int*)b - *(long int*)a); // Sort in descending order
}

// Comparison function for qsort (ascending order)
int compare_longs(const void *a, const void *b) {
    return (*(long int *)a < *(long int *)b) ? -1 : (*(long int *)a > *(long int *)b);
}


// smooths counts if smooth parameter is set to 2
long int smooth_kmer_counts(long int *******results,
                          SingleAndTripleGapKmer *g,
                          short int file_number,
                          short int smooth) {
    
    // Return 0 if memory is not allocated
    if (Kmer_out_of_range_linear(results, g, file_number)) return 0;
    
    long int kmer = g->triplegap_kmer_value;
    
    // Return unsmoothed count if smooth is 0
    if (smooth != 2) {
        /* int dummy;   // debugging uninitialized value from stack allocation
        if (file_number > -999) dummy = 1;
        if (g->original_kmer_number > -999) dummy = 1;
        if (g->added_length > -999) dummy = 1;
        if (g->gap_pos[0] > -999) dummy = 1;
        if (g->gap_pos[1] > -999) dummy = 1;
        if (g->gap_pos[2] > -999) dummy = 1;
        if (kmer > -999) dummy = 1; */
        return results[file_number][g->original_kmer_number][g->added_length][g->gap_pos[0]][g->gap_pos[1]][g->gap_pos[2]][kmer];
    }
    
    // Get original kmer count
    long int original_count = results[file_number][g->original_kmer_number][g->added_length][g->gap_pos[0]][g->gap_pos[1]][g->gap_pos[2]][kmer];
    
    // Store neighbor counts (excluding original kmer)
    long int neighbor_counts[3 * g->added_length]; // Max neighbors: 3 options per position
    int count = 0;
    
    // Generate all Hamming distance 1 neighbors
    long int lowbit = 1;  // For A-C or G-T transitions
    long int highbit = 2; // For A-G or C-T transitions
    long int compared_kmer;
    
    for(short int position = 0; position < g->added_length; position++, lowbit <<= 2, highbit <<= 2) {
        // First substitution (A↔C, G↔T)
        compared_kmer = lowbit ^ kmer;
        neighbor_counts[count++] = results[file_number][g->original_kmer_number][g->added_length][g->gap_pos[0]][g->gap_pos[1]][g->gap_pos[2]][compared_kmer];
        
        // Second substitution (A↔G, C↔T)
        compared_kmer = highbit ^ kmer;
        neighbor_counts[count++] = results[file_number][g->original_kmer_number][g->added_length][g->gap_pos[0]][g->gap_pos[1]][g->gap_pos[2]][compared_kmer];
        
        // Third substitution (A↔T, C↔G)
        compared_kmer = lowbit ^ compared_kmer;
        neighbor_counts[count++] = results[file_number][g->original_kmer_number][g->added_length][g->gap_pos[0]][g->gap_pos[1]][g->gap_pos[2]][compared_kmer];
    }
    
    // Sort counts in descending order
    qsort(neighbor_counts, count, sizeof(long int), compare_counts);
    
    // Average top 10 neighbor counts (or all if less than 10)
    int top_n = (count < 10) ? count : 10;
    long int sum = 0;
    for(int i = 0; i < top_n; i++) {
        sum += neighbor_counts[i];
    }
    long int neighbor_average = sum / top_n;
    
    // Return 50-50 blend of original count and neighbor average
    return (original_count + neighbor_average) / 2;
}

// inline function to check if another kmer count is higher than current kmer
static inline int check_kmer(
    long int *******results,
    SingleAndTripleGapKmer *orig_kmer,
    short int file_number,
    short int current_kmer_length,
    short int current_gap_position,
    short int current_gap_length,
    long int compared_kmer,
    long int kmer1_incidence,
    int background_subtraction,
    double scale_factor,
    double kmer_length_difference_cutoff,
    short int original_kmer_length,
    short int robust)  // Added parameter
{
    
    MergedKmer full_kmer;
    full_kmer.kmer_value = compared_kmer;
    full_kmer.kmer_length = current_kmer_length;
    full_kmer.gap_position = current_gap_position;
    full_kmer.gap_length = current_gap_length;
    
    if (orig_kmer->debug > 2) {
        printf(" ");
        Kmerprint(compared_kmer, current_kmer_length, current_gap_position, current_gap_length);
        printf(" ");
        Kmerprint(orig_kmer->original_kmer_value, orig_kmer->original_kmer_length, orig_kmer->original_gap_position, orig_kmer->original_gap_length);
        printf(" ");
        Triplegap_kmerprint(orig_kmer);
        printf(" Count: %li", (long int) kmer1_incidence);
        
    }
    
    // convert index from 5d to 7d
    short int map_index = Align_to_identify_triplegapkmer(orig_kmer, &full_kmer);
    if (map_index == 0) return 0;
    
    long int kmer2_incidence = smooth_kmer_counts(results, orig_kmer, file_number, background_subtraction);
    if (robust) kmer2_incidence += (int) sqrt(kmer2_incidence);
    
    long int background;

    if (orig_kmer->debug == 1) {
        printf(" **");
        Kmerprint(full_kmer.kmer_value, full_kmer.kmer_length, full_kmer.gap_position, full_kmer.gap_length);
        printf(" ");
        Kmerprint(orig_kmer->original_kmer_value, orig_kmer->original_kmer_length, orig_kmer->original_gap_position, orig_kmer->original_gap_length);
        printf(" pos:len %i:%i,%i:%i,%i:%i ", orig_kmer->gap_pos[0], orig_kmer->gap_len[0], orig_kmer->gap_pos[1], orig_kmer->gap_len[1], orig_kmer->gap_pos[2], orig_kmer->gap_len[2]);
        Triplegap_kmerprint(orig_kmer);
        printf(" Counts: %li vs %li, index %i **\t", (long int) kmer1_incidence, (long int) kmer2_incidence, map_index);
        
    }
    

    
    if (background_subtraction) {
        if (scale_factor >= 1.0) {
            background = smooth_kmer_counts(results, orig_kmer, 0, background_subtraction);
            kmer2_incidence -= background * (long int)scale_factor;
        } else {
            background = smooth_kmer_counts(results, orig_kmer, 0, background_subtraction);
            kmer2_incidence = kmer2_incidence * (long int)(1.0/scale_factor) - background;
        }
    }
    
    // Regular comparison for same length kmers
    if (current_kmer_length == original_kmer_length) {
        if (orig_kmer->debug == -1 && kmer2_incidence == kmer1_incidence && kmer1_incidence > 1000) {
            printf (" self? ");
            Kmerprint(full_kmer.kmer_value, full_kmer.kmer_length, full_kmer.gap_position, full_kmer.gap_length);
            printf(" ");
            Triplegap_kmerprint(orig_kmer);
        }
        return (kmer2_incidence > kmer1_incidence);
    }
    // For longer kmers, divide kmer1 by cutoff
    else if (current_kmer_length > original_kmer_length) {
        return (kmer2_incidence > kmer1_incidence * kmer_length_difference_cutoff);
    }
    // For shorter kmers, multiply kmer2 by cutoff
    else {
        return (kmer2_incidence * kmer_length_difference_cutoff > kmer1_incidence);
    }
}

/* SUBROUTINE THAT DETERMINES IF GAP IS AT EITHER OF THE CENTER POSITIONS, IF count_also is set to != 1 returns true */
short int Centergap (short int count_also_spaced_kmers, short int kmer_length, short int gap_position)
{
if (count_also_spaced_kmers != 1) return (1);
if (gap_position == kmer_length / 2) return (1);
if (kmer_length % 2 == 1 && gap_position == kmer_length / 2 - 1) return (1);
else return (0);
}

/* Localmax algorithm adapted to the merged kmer approach */
short int Mergedkmer_Localmax(long int *******results, SingleAndTripleGapKmer *gap_iter, MergedKmer *merged_kmer,
                            short int background_subtraction,
                            double scale_factor, short int robust, double kmer_length_difference_cutoff) {

    // only consider signal counts
    short int file_number = 1;
    
   long int current_kmer = merged_kmer->kmer_value;
   short int current_kmer_length = merged_kmer->kmer_length;
   short int shortest_kmer = 1;
    short int longest_kmer_counted = gap_iter->original_kmer_length + gap_iter->max_length; // input here are the extended kmers, so original gapped kmer + max added length input parameter is the correct value;
   short int current_gap_position = merged_kmer->gap_position;
   short int current_gap_length = merged_kmer->gap_length;
    
    short int original_kmer_length = current_kmer_length;
    short int count_also_spaced_kmers = 2;
    short int too_long_kmer = longest_kmer_counted + 1;

    // Calculate initial kmer1_incidence
    long int kmer1_incidence;
    if (background_subtraction) {
        long int background = smooth_kmer_counts(results, gap_iter, 0, background_subtraction);
        if (scale_factor >= 1.0) {
            kmer1_incidence = smooth_kmer_counts(results, gap_iter, file_number, background_subtraction)
                             - background * (long int)scale_factor;
        } else {
            kmer1_incidence = smooth_kmer_counts(results, gap_iter, file_number, background_subtraction)
                             * (long int)(1.0/scale_factor) - background;
        }
    } else {
        kmer1_incidence = smooth_kmer_counts(results, gap_iter, file_number, background_subtraction);
    }

    //if (robust) kmer1_incidence -= sqrt(kmer1_incidence);
    
    long int compared_kmer = current_kmer;
    signed short int position;
    short int counter;
    short int first_half;
    long int lowbit = 1;
    long int highbit = 2;
    short int shift;
    short int true_gap_position = current_gap_position;
    short int true_gap_length = current_gap_length;
    short int start;
    short int end;
    short int left = 0;
    short int right = 1;
    signed short int position2;

    short int number_of_comparisons;

    //gap_iter->debug = 0;
    if (gap_iter->debug == 1) printf (" len:%i", current_kmer_length);
    /* Substitution; HAMMING OF 1 */
    for(position=0, number_of_comparisons = 0; position < current_kmer_length; position++, lowbit <<= 2, highbit <<= 2) {
        
        compared_kmer = lowbit ^ current_kmer;
        if(check_kmer(results, gap_iter, file_number, current_kmer_length, current_gap_position, current_gap_length, compared_kmer, kmer1_incidence, background_subtraction, scale_factor, kmer_length_difference_cutoff, original_kmer_length, robust)) return 0;
        
        compared_kmer = highbit ^ current_kmer;
        if(check_kmer(results, gap_iter, file_number, current_kmer_length, current_gap_position, current_gap_length, compared_kmer, kmer1_incidence, background_subtraction, scale_factor, kmer_length_difference_cutoff, original_kmer_length, robust)) return 0;
        
        //gap_iter->debug = -1;
        compared_kmer = lowbit ^ compared_kmer;
        if(check_kmer(results, gap_iter, file_number, current_kmer_length, current_gap_position, current_gap_length, compared_kmer, kmer1_incidence, background_subtraction, scale_factor, kmer_length_difference_cutoff, original_kmer_length, robust)) return 0;
        number_of_comparisons++;

    }
    //gap_iter->debug = 0;
    
    /* Shift */
    shift = (current_gap_length != 0);
    current_gap_position += shift;
    
    if (current_gap_length == 0 || (count_also_spaced_kmers == 2 - SHIFTED_GAP_POSITIONS && current_gap_position < current_kmer_length)
        || (current_kmer_length % 2 == 1 && count_also_spaced_kmers == 1 && current_gap_position == current_kmer_length / 2 + 1)) {
        
        compared_kmer = (current_kmer >> 2) & lowmask_ULL[current_kmer_length-1];
        if(check_kmer(results, gap_iter, file_number, current_kmer_length, current_gap_position, current_gap_length, compared_kmer, kmer1_incidence, background_subtraction, scale_factor, kmer_length_difference_cutoff, original_kmer_length, robust)) return 0;
        
        lowbit >>= 2; highbit >>= 2;
        compared_kmer = lowbit ^ compared_kmer;
        if(check_kmer(results, gap_iter, file_number, current_kmer_length, current_gap_position, current_gap_length, compared_kmer, kmer1_incidence, background_subtraction, scale_factor, kmer_length_difference_cutoff, original_kmer_length, robust)) return 0;
        
        compared_kmer = highbit ^ compared_kmer;
        if(check_kmer(results, gap_iter, file_number, current_kmer_length, current_gap_position, current_gap_length, compared_kmer, kmer1_incidence, background_subtraction, scale_factor, kmer_length_difference_cutoff, original_kmer_length, robust)) return 0;
        
        compared_kmer = lowbit ^ compared_kmer;
        if(check_kmer(results, gap_iter, file_number, current_kmer_length, current_gap_position, current_gap_length, compared_kmer, kmer1_incidence, background_subtraction, scale_factor, kmer_length_difference_cutoff, original_kmer_length, robust)) return 0;
    }
    
    current_gap_position -= shift;
    current_gap_position -= shift;
    if (current_gap_length == 0 || (count_also_spaced_kmers == 2 - SHIFTED_GAP_POSITIONS && current_gap_position > 0)
        || (current_kmer_length % 2 == 1 && count_also_spaced_kmers == 1 && current_gap_position == current_kmer_length / 2)) {
        
        compared_kmer = (current_kmer << 2) & lowmask_ULL[current_kmer_length-1];
        if(check_kmer(results, gap_iter, file_number, current_kmer_length, current_gap_position, current_gap_length, compared_kmer, kmer1_incidence, background_subtraction, scale_factor, kmer_length_difference_cutoff, original_kmer_length, robust)) return 0;
        
        lowbit = 1; highbit = 2;
        compared_kmer = lowbit ^ compared_kmer;
        if(check_kmer(results, gap_iter, file_number, current_kmer_length, current_gap_position, current_gap_length, compared_kmer, kmer1_incidence, background_subtraction, scale_factor, kmer_length_difference_cutoff, original_kmer_length, robust)) return 0;
        
        compared_kmer = highbit ^ compared_kmer;
        if(check_kmer(results, gap_iter, file_number, current_kmer_length, current_gap_position, current_gap_length, compared_kmer, kmer1_incidence, background_subtraction, scale_factor, kmer_length_difference_cutoff, original_kmer_length, robust)) return 0;
        
        compared_kmer = lowbit ^ compared_kmer;
        if(check_kmer(results, gap_iter, file_number, current_kmer_length, current_gap_position, current_gap_length, compared_kmer, kmer1_incidence, background_subtraction, scale_factor, kmer_length_difference_cutoff, original_kmer_length, robust)) return 0;
    }
    
    current_gap_position += shift;
    lowbit = 1; highbit = 2;
    if (count_also_spaced_kmers != 0) {
        if(current_gap_position == 0) current_gap_position = current_kmer_length / 2;
        
        /* Longer Gap */
        current_gap_length++;
        if (current_gap_length < Nlength - current_kmer_length && current_gap_length < MAX_GAP_LENGTH && true_gap_length != 0) {
            compared_kmer = (current_kmer & highmask_ULL[current_kmer_length-current_gap_position]) | ((current_kmer << 2) & (lowmask_ULL[current_kmer_length-current_gap_position-1]));
            
            for(first_half = 0; first_half < 2; first_half++) {
                if(check_kmer(results, gap_iter, file_number, current_kmer_length, current_gap_position, current_gap_length, compared_kmer, kmer1_incidence, background_subtraction, scale_factor, kmer_length_difference_cutoff, original_kmer_length, robust)) return 0;
                
                compared_kmer = lowbit ^ compared_kmer;
                if(check_kmer(results, gap_iter, file_number, current_kmer_length, current_gap_position, current_gap_length, compared_kmer, kmer1_incidence, background_subtraction, scale_factor, kmer_length_difference_cutoff, original_kmer_length, robust)) return 0;
                
                compared_kmer = highbit ^ compared_kmer;
                if(check_kmer(results, gap_iter, file_number, current_kmer_length, current_gap_position, current_gap_length, compared_kmer, kmer1_incidence, background_subtraction, scale_factor, kmer_length_difference_cutoff, original_kmer_length, robust)) return 0;
                
                compared_kmer = lowbit ^ compared_kmer;
                if(check_kmer(results, gap_iter, file_number, current_kmer_length, current_gap_position, current_gap_length, compared_kmer, kmer1_incidence, background_subtraction, scale_factor, kmer_length_difference_cutoff, original_kmer_length, robust)) return 0;
                
                compared_kmer = ((current_kmer >> 2) & highmask_ULL[current_kmer_length-current_gap_position]) | ((current_kmer) & (lowmask_ULL[current_kmer_length-current_gap_position-1]));
                lowbit <<= ((current_kmer_length - 1)*2); highbit <<= ((current_kmer_length - 1)*2);
            }
        }
        
        /* Shorter Gap */
        current_gap_length--;
        current_gap_length--;
        
        lowbit = 1; highbit = 2;
        lowbit <<= ((current_kmer_length - true_gap_position - 1)*2); highbit <<= ((current_kmer_length - true_gap_position - 1)*2);
        if (current_gap_length == 0) current_gap_position = 0;
        if (current_gap_length >= 0) {
            compared_kmer = (current_kmer & highmask_ULL[current_kmer_length - true_gap_position]) | ((current_kmer >> 2) & (lowmask_ULL[current_kmer_length-true_gap_position-1]));
            
            for(first_half = 0; first_half < 2; first_half++) {
                if(check_kmer(results, gap_iter, file_number, current_kmer_length, current_gap_position, current_gap_length, compared_kmer, kmer1_incidence, background_subtraction, scale_factor, kmer_length_difference_cutoff, original_kmer_length, robust)) return 0;
                
                compared_kmer = lowbit ^ compared_kmer;
                if(check_kmer(results, gap_iter, file_number, current_kmer_length, current_gap_position, current_gap_length, compared_kmer, kmer1_incidence, background_subtraction, scale_factor, kmer_length_difference_cutoff, original_kmer_length, robust)) return 0;
                
                compared_kmer = highbit ^ compared_kmer;
                if(check_kmer(results, gap_iter, file_number, current_kmer_length, current_gap_position, current_gap_length, compared_kmer, kmer1_incidence, background_subtraction, scale_factor, kmer_length_difference_cutoff, original_kmer_length, robust)) return 0;
                
                compared_kmer = lowbit ^ compared_kmer;
                if(check_kmer(results, gap_iter, file_number, current_kmer_length, current_gap_position, current_gap_length, compared_kmer, kmer1_incidence, background_subtraction, scale_factor, kmer_length_difference_cutoff, original_kmer_length, robust)) return 0;
                
                compared_kmer = lowmask_ULL[current_kmer_length-1] & ((current_kmer << 2) & highmask_ULL[current_kmer_length - true_gap_position]) | ((current_kmer) & (lowmask_ULL[current_kmer_length-true_gap_position-1]));
                lowbit <<= 2; highbit <<= 2;
            }
        }
        /* Different gap positions */
                current_gap_length = true_gap_length;
                if ((count_also_spaced_kmers == 2 || (count_also_spaced_kmers == 1 && current_kmer_length % 2 == 1)) && true_gap_length > 0) {
                    current_gap_position = true_gap_position + 1;
                    lowbit = 1; highbit = 2;
                    lowbit <<= ((current_kmer_length - current_gap_position)*2); highbit <<= ((current_kmer_length - current_gap_position)*2);
                    compared_kmer = (current_kmer & highmask_ULL[current_kmer_length - current_gap_position]) | ((current_kmer) & (lowmask_ULL[current_kmer_length-current_gap_position-1]));

                    for(first_half = 0; first_half < 2; first_half++) {
                        if (current_gap_position < current_kmer_length && current_gap_position > 0 && (count_also_spaced_kmers == 2 ||
                            (count_also_spaced_kmers == 1 && ((current_gap_position == current_kmer_length / 2) || current_gap_position == current_kmer_length / 2 + 1)))) {
                            
                            if(check_kmer(results, gap_iter, file_number, current_kmer_length, current_gap_position, current_gap_length, compared_kmer, kmer1_incidence, background_subtraction, scale_factor, kmer_length_difference_cutoff, original_kmer_length, robust)) return 0;
                            compared_kmer = lowbit ^ compared_kmer;
                            if(check_kmer(results, gap_iter, file_number, current_kmer_length, current_gap_position, current_gap_length, compared_kmer, kmer1_incidence, background_subtraction, scale_factor, kmer_length_difference_cutoff, original_kmer_length, robust)) return 0;
                            compared_kmer = highbit ^ compared_kmer;
                            if(check_kmer(results, gap_iter, file_number, current_kmer_length, current_gap_position, current_gap_length, compared_kmer, kmer1_incidence, background_subtraction, scale_factor, kmer_length_difference_cutoff, original_kmer_length, robust)) return 0;
                            compared_kmer = lowbit ^ compared_kmer;
                            if(check_kmer(results, gap_iter, file_number, current_kmer_length, current_gap_position, current_gap_length, compared_kmer, kmer1_incidence, background_subtraction, scale_factor, kmer_length_difference_cutoff, original_kmer_length, robust)) return 0;
                        }
                        
                        current_gap_position--;
                        current_gap_position--;
                        compared_kmer = ((current_kmer) & highmask_ULL[current_kmer_length - current_gap_position]) | ((current_kmer) & (lowmask_ULL[current_kmer_length-current_gap_position-1]));
                        lowbit <<= 2; highbit <<= 2;
                    }
                }
                
                start = 1;
                end = current_kmer_length;
                
                /* Compare ungapped kmer to all single gaps */
                if (count_also_spaced_kmers != 0 && true_gap_length == 0) {
                    if(count_also_spaced_kmers == 1) {
                        start = current_kmer_length / 2;
                        end = start + 1 + (current_kmer_length % 2);
                    }
                    
                    for(current_gap_position = start, current_gap_length = 1; current_gap_position < end; current_gap_position++) {
                        compared_kmer = (current_kmer & highmask_ULL[current_kmer_length - current_gap_position]) | ((current_kmer << 2) & (lowmask_ULL[current_kmer_length-current_gap_position-1]));
                        
                        for(lowbit = 1, highbit = 2, first_half = 0; first_half < 2; first_half++) {
                            if(check_kmer(results, gap_iter, file_number, current_kmer_length, current_gap_position, current_gap_length, compared_kmer, kmer1_incidence, background_subtraction, scale_factor, kmer_length_difference_cutoff, original_kmer_length, robust)) return 0;
                            compared_kmer = lowbit ^ compared_kmer;
                            if(check_kmer(results, gap_iter, file_number, current_kmer_length, current_gap_position, current_gap_length, compared_kmer, kmer1_incidence, background_subtraction, scale_factor, kmer_length_difference_cutoff, original_kmer_length, robust)) return 0;
                            compared_kmer = highbit ^ compared_kmer;
                            if(check_kmer(results, gap_iter, file_number, current_kmer_length, current_gap_position, current_gap_length, compared_kmer, kmer1_incidence, background_subtraction, scale_factor, kmer_length_difference_cutoff, original_kmer_length, robust)) return 0;
                            compared_kmer = lowbit ^ compared_kmer;
                            if(check_kmer(results, gap_iter, file_number, current_kmer_length, current_gap_position, current_gap_length, compared_kmer, kmer1_incidence, background_subtraction, scale_factor, kmer_length_difference_cutoff, original_kmer_length, robust)) return 0;
                            
                            compared_kmer = ((current_kmer >> 2) & highmask_ULL[current_kmer_length - current_gap_position]) | ((current_kmer) & (lowmask_ULL[current_kmer_length-current_gap_position-1]));
                            lowbit <<= ((current_kmer_length - 1)*2); highbit <<= ((current_kmer_length - 1)*2);
                        }
                    }
                }
            }

            /* Compare with shorter kmer */
            current_gap_position = true_gap_position;
            current_gap_length = true_gap_length;
            current_kmer_length--;
            end = current_kmer_length;
            if (current_kmer_length >= shortest_kmer) {
                if(current_gap_length == 0) {
                    if (count_also_spaced_kmers != 0) {
                        if(count_also_spaced_kmers == 1) {
                            start = current_kmer_length / 2;
                            end = start + 1 + (current_kmer_length % 2);
                        }
                        for(current_gap_position = start, current_gap_length = 1; current_gap_position < end; current_gap_position++) {
                            compared_kmer = ((current_kmer >> 2) & highmask_ULL[current_kmer_length - current_gap_position]) | ((current_kmer) & (lowmask_ULL[current_kmer_length-current_gap_position-1]));
                            if(check_kmer(results, gap_iter, file_number, current_kmer_length, current_gap_position, current_gap_length, compared_kmer, kmer1_incidence, background_subtraction, scale_factor, kmer_length_difference_cutoff, original_kmer_length, robust)) return 0;
                        }
                    }
                }

                current_gap_position = true_gap_position;
                current_gap_length = true_gap_length;
                
                if (count_also_spaced_kmers != 1) {left = 1; right = 1;}
                if (current_gap_position == current_kmer_length / 2 && current_kmer_length % 2 == 0 && count_also_spaced_kmers == 1) {left = 1; right = 0;}
                if (current_kmer_length % 2 == 1 && count_also_spaced_kmers == 1) left = 1;
                
                /* Left part */
                if (current_gap_position < current_kmer_length) {
                    if (left == 1 || true_gap_length == 0) {
                        compared_kmer = (current_kmer >> 2);
                        if(check_kmer(results, gap_iter, file_number, current_kmer_length, current_gap_position, current_gap_length, compared_kmer, kmer1_incidence, background_subtraction, scale_factor, kmer_length_difference_cutoff, original_kmer_length, robust)) return 0;
                    }
                }
                
                /* Right part */
                if (current_gap_position != 1 && right == 1) {
                    if (current_gap_position > 0) current_gap_position--;
                    compared_kmer = (current_kmer & lowmask_ULL[current_kmer_length-1]);
                    if(check_kmer(results, gap_iter, file_number, current_kmer_length, current_gap_position, current_gap_length, compared_kmer, kmer1_incidence, background_subtraction, scale_factor, kmer_length_difference_cutoff, original_kmer_length, robust)) return 0;
                }
                
                current_gap_position = true_gap_position;
                
                /* Shorter with longer gap */
                if(count_also_spaced_kmers != 0 && true_gap_position != 0 && current_gap_length < MAX_GAP_LENGTH) {
                    current_gap_length++;
                    if (current_gap_position < current_kmer_length && left == 1) {
                        compared_kmer = ((current_kmer >> 2) & highmask_ULL[current_kmer_length - current_gap_position]) | ((current_kmer) & (lowmask_ULL[current_kmer_length-current_gap_position-1]));
                        if(check_kmer(results, gap_iter, file_number, current_kmer_length, current_gap_position, current_gap_length, compared_kmer, kmer1_incidence, background_subtraction, scale_factor, kmer_length_difference_cutoff, original_kmer_length, robust)) return 0;
                    }
                    
                    if(current_gap_position > 1 && right == 1) {
                        current_gap_position--;
                        compared_kmer = ((current_kmer >> 2) & highmask_ULL[current_kmer_length - current_gap_position]) | ((current_kmer) & (lowmask_ULL[current_kmer_length-current_gap_position-1]));
                        if(check_kmer(results, gap_iter, file_number, current_kmer_length, current_gap_position, current_gap_length, compared_kmer, kmer1_incidence, background_subtraction, scale_factor, kmer_length_difference_cutoff, original_kmer_length, robust)) return 0;
                    }
                    current_gap_length--;
                    current_gap_position++;
                }

                /* Compare hanging single base to ungapped kmer */
                current_gap_position = true_gap_position;
                current_gap_length = true_gap_length;
                if(current_gap_position == 1) {
                    current_gap_length = 0;
                    current_gap_position = 0;
                    compared_kmer = (current_kmer & lowmask_ULL[current_kmer_length-1]);
                    if(check_kmer(results, gap_iter, file_number, current_kmer_length, current_gap_position, current_gap_length, compared_kmer, kmer1_incidence, background_subtraction, scale_factor, kmer_length_difference_cutoff, original_kmer_length, robust)) return 0;
                }
                else if(current_kmer_length-current_gap_position == 0) {
                    current_gap_length = 0;
                    current_gap_position = 0;
                    compared_kmer = (current_kmer >> 2);
                    if(check_kmer(results, gap_iter, file_number, current_kmer_length, current_gap_position, current_gap_length, compared_kmer, kmer1_incidence, background_subtraction, scale_factor, kmer_length_difference_cutoff, original_kmer_length, robust)) return 0;
                }
            }

            current_gap_position = true_gap_position;
            current_gap_length = true_gap_length;
            
            /* Compare with longer kmer */
            current_kmer_length++;
            current_kmer_length++;
            if (current_kmer_length < too_long_kmer) {
                compared_kmer = (current_kmer << 2);
                for(lowbit = 1, highbit = 2, first_half = 0; first_half < 2; first_half++) {
                    if(count_also_spaced_kmers != 1 || true_gap_length == 0 ||
                       (first_half == 1 || (count_also_spaced_kmers == 1 && current_gap_position == current_kmer_length / 2))) {
                        
                        //gap_iter->debug = 1;
                        if(check_kmer(results, gap_iter, file_number, current_kmer_length, current_gap_position, current_gap_length, compared_kmer, kmer1_incidence, background_subtraction, scale_factor, kmer_length_difference_cutoff, original_kmer_length, robust)) return 0;
                        compared_kmer = lowbit ^ compared_kmer;
                        if(check_kmer(results, gap_iter, file_number, current_kmer_length, current_gap_position, current_gap_length, compared_kmer, kmer1_incidence, background_subtraction, scale_factor, kmer_length_difference_cutoff, original_kmer_length, robust)) return 0;
                        compared_kmer = highbit ^ compared_kmer;
                        if(check_kmer(results, gap_iter, file_number, current_kmer_length, current_gap_position, current_gap_length, compared_kmer, kmer1_incidence, background_subtraction, scale_factor, kmer_length_difference_cutoff, original_kmer_length, robust)) return 0;
                        compared_kmer = lowbit ^ compared_kmer;
                        if(check_kmer(results, gap_iter, file_number, current_kmer_length, current_gap_position, current_gap_length, compared_kmer, kmer1_incidence, background_subtraction, scale_factor, kmer_length_difference_cutoff, original_kmer_length, robust)) return 0;
                        //gap_iter->debug = 0;
                    }
                    
                    compared_kmer = current_kmer;
                    lowbit <<= ((current_kmer_length - 1)*2); highbit <<= ((current_kmer_length - 1)*2);
                    if (true_gap_length == 0) continue;
                    current_gap_position++;
                    if (current_gap_position > current_kmer_length || (count_also_spaced_kmers == 1 && current_gap_position != (current_kmer_length / 2 + (current_kmer_length % 2)))) break;
                }
                
                current_gap_position = true_gap_position;
                current_gap_length = true_gap_length;
                /* Longer with shorter gap */
                //gap_iter->debug = 1;
                        if(count_also_spaced_kmers != 0 && true_gap_length >= 1) {
                            current_gap_length--;
                            lowbit = 1; highbit = 2;
                            lowbit <<= ((current_kmer_length - current_gap_position-1)*2); highbit <<= ((current_kmer_length - current_gap_position-1)*2);
                            compared_kmer = ((current_kmer << 2) & highmask_ULL[current_kmer_length - current_gap_position]) | ((current_kmer) & (lowmask_ULL[current_kmer_length-current_gap_position-1]));
                            if (current_gap_length == 0) current_gap_position = 0;
                            
                            for(first_half = 0; first_half < 2; first_half++) {
                                if (current_gap_position < current_kmer_length && (current_gap_position == 0 || Centergap (count_also_spaced_kmers, current_kmer_length, current_gap_position))) {
                                    if(check_kmer(results, gap_iter, file_number, current_kmer_length, current_gap_position, current_gap_length, compared_kmer, kmer1_incidence, background_subtraction, scale_factor, kmer_length_difference_cutoff, original_kmer_length, robust)) return 0;
                                    compared_kmer = lowbit ^ compared_kmer;
                                    if(check_kmer(results, gap_iter, file_number, current_kmer_length, current_gap_position, current_gap_length, compared_kmer, kmer1_incidence, background_subtraction, scale_factor, kmer_length_difference_cutoff, original_kmer_length, robust)) return 0;
                                    compared_kmer = highbit ^ compared_kmer;
                                    if(check_kmer(results, gap_iter, file_number, current_kmer_length, current_gap_position, current_gap_length, compared_kmer, kmer1_incidence, background_subtraction, scale_factor, kmer_length_difference_cutoff, original_kmer_length, robust)) return 0;
                                    compared_kmer = lowbit ^ compared_kmer;
                                    if(check_kmer(results, gap_iter, file_number, current_kmer_length, current_gap_position, current_gap_length, compared_kmer, kmer1_incidence, background_subtraction, scale_factor, kmer_length_difference_cutoff, original_kmer_length, robust)) return 0;
                                }
                                
                                current_gap_position++;
                                if (current_gap_position > current_kmer_length || current_gap_position == 1) break;
                                compared_kmer = ((current_kmer << 2) & highmask_ULL[current_kmer_length - current_gap_position]) | ((current_kmer) & (lowmask_ULL[current_kmer_length-current_gap_position-1]));
                            }
                            current_gap_length++;
                        }
                    //gap_iter->debug = 0;
                    }
                    current_gap_position = true_gap_position;
                    return 1;
                }

void Test_Compute_triplegap_kmer(SingleAndTripleGapKmer* o, MergedKmer* merged_kmer) {
Identify_triplegapkmer(o, merged_kmer);
}

// iterator that finds allowed gap patterns for a triplegapped fill kmer straddling a gapped kmer
short int move_gaps(SingleAndTripleGapKmer* iter) {

    // reiterates until allowed pattern found
    for (;;)
    {
        
        if((*iter).gap_pos[2] < (*iter).added_length) { ((*iter).gap_pos[2])++;}
        
        else if ((*iter).gap_pos[1] < (*iter).added_length) { ((*iter).gap_pos[1])++; (*iter).gap_pos[2] = (*iter).gap_pos[1];}
        
        else if ((*iter).gap_pos[0] < (*iter).added_length) { ((*iter).gap_pos[0])++; (*iter).gap_pos[1] = (*iter).gap_pos[0]; (*iter).gap_pos[2] = (*iter).gap_pos[0];}
        
        else if ((*iter).added_length >= (*iter).max_length) return 0;
        
        else {
            (*iter).added_length++;
            for (int i = 0; i < 3; i++) (*iter).gap_pos[i] = 0;
            }
        
        // allowed patterns have fewer or equal number of bases within the original gap. All cases with fewer bases are allowed.
        if ((*iter).gap_pos[2] - (*iter).gap_pos[0] < (*iter).original_gap_length) break;
        // but only one case where the gap is completely filled is allowed
        if ((*iter).gap_pos[2] - (*iter).gap_pos[0] == (*iter).original_gap_length && (*iter).gap_pos[1] == (*iter).gap_pos[0]) break;
        
    }
    
    //
    //Triplegap_kmerprint(iter);
    
    // calculates filler triplegap kmer lengths based on original kmer
    // first gap is always as long as the first half of the original kmer
    (*iter).gap_len[0] = (*iter).original_gap_position;
    
    // second gap is shorter than original gap length (sum of lengths of the filler kmer bases in the gap)
    (*iter).gap_len[1] = (*iter).original_gap_length + (*iter).gap_pos[0] - (*iter).gap_pos[2];
    
    // last gap must end at original kmer lenght + original gap length (total original kmer width)
    // (*iter).gap_len[2] = (*iter).original_kmer_length + (*iter).original_gap_length - (*iter).gap_len[0] - (*iter).gap_len[1] + (*iter).gap_pos[0] - (*iter).gap_pos[2];
    
    // and more simply, last gap jumps over the right part of the original kmer
    (*iter).gap_len[2] = (*iter).original_kmer_length - (*iter).original_gap_position;
    // (*iter).gap_len[2] = (*iter).original_kmer_length + (*iter).original_gap_length - (*iter).gap_len[0] - (*iter).gap_len[1] + (*iter).gap_pos[0] - (*iter).gap_pos[2];
    
    
    //printf("\nOrig kmer ");
    //Kmerprint(iter->original_kmer_value, iter->original_kmer_length, iter->original_gap_position, iter->original_gap_length); fflush(stdout);
    //printf(":added lengths-%i:%i gap pos:lengths %i:%i,%i:%i,%i:%i", (*iter).original_kmer_length, (*iter).added_length, (*iter).gap_pos[0], (*iter).gap_len[0], (*iter).gap_pos[1], (*iter).gap_len[1], (*iter).gap_pos[2], (*iter).gap_len[2]);
    
    return 1;
    }

void Set_original_kmer_value(const char* sequence, SingleAndTripleGapKmer* gap_info) {
   gap_info->original_kmer_value = 0;
   
   for(int i = 0; i < strlen(sequence); i++) {
       if(sequence[i] != 'n') {
           gap_info->original_kmer_value = gap_info->original_kmer_value << 2;
           switch(sequence[i]) {
               case 'A': gap_info->original_kmer_value |= 0; break;
               case 'C': gap_info->original_kmer_value |= 1; break;
               case 'G': gap_info->original_kmer_value |= 2; break;
               case 'T': gap_info->original_kmer_value |= 3; break;
           }
       }
   }
}

// Get gap info from sequence and set kmer value
void Set_original_gap_info_and_init(const char* sequence, SingleAndTripleGapKmer* g, short int max_length, short int input_kmer_number) {
    g->added_length = 0;
    g->max_length = max_length;
    g->original_gap_position = 0;
    g->original_gap_length = 0;
    g->original_kmer_value = 0;
    g->original_kmer_number = input_kmer_number;
    g->debug = 0;
    
    for(int i = 0; i < 3; i++) {
        g->gap_pos[i] = 0;
        g->gap_len[i] = 0;
    }
    
    // First count length and find gap
    int seq_len = strlen(sequence);
    for(int i = 0; i < seq_len; i++) {
        if(sequence[i] == 'n') {
            if(g->original_gap_length == 0) {
                g->original_gap_position = i;
            }
            (g->original_gap_length)++;
        }
    }
    g->original_kmer_length = seq_len - g->original_gap_length;
    
    Set_original_kmer_value(sequence, g);
    // Build value using bit shifts

    // Kmerprint(g->original_kmer_value, g->original_kmer_length, g->original_gap_position, g->original_gap_length); fflush(stdout);
}

long int******* allocate_counts(int num_files, LocalMaxTable* table, int max_extension) {
    
    int num_local_max = table->num_entries;
    long total_bytes = 0;
    
    // Dimension 1: top-level array of pointers
    long int******* counts = calloc(num_files, sizeof(*counts));
    if (!counts) {
        fprintf(stderr, "Failed to allocate counts\n");
        exit(EXIT_FAILURE);
    }
    total_bytes += num_files * sizeof(*counts);

    for (int f = 0; f < num_files; f++) {
        // Dimension 2
        if (!counts[f]) {
            counts[f] = calloc(num_local_max, sizeof(*counts[f]));
            if (!counts[f]) {
                fprintf(stderr, "Failed to allocate counts[%d]\n", f);
                exit(EXIT_FAILURE);
            }
            total_bytes += num_local_max * sizeof(*counts[f]);
        }

        for (int l = 0; l < num_local_max; l++) {
            // Dimension 3
            if (!counts[f][l]) {
                counts[f][l] = calloc(MAX_ADDED_LEN + 1, sizeof(*counts[f][l]));
                if (!counts[f][l]) {
                    fprintf(stderr, "Failed to allocate counts[%d][%d]\n", f, l);
                    exit(EXIT_FAILURE);
                }
                total_bytes += (MAX_ADDED_LEN + 1) * sizeof(*counts[f][l]);
            }

            // Initialize gap iterator for this local max entry
            SingleAndTripleGapKmer gap_iter;
            Set_original_gap_info_and_init(table->entries[l].sequence, &gap_iter, max_extension, l);

            // printf("\nAllocating for local max: %s\n", table->entries[l].sequence);
            // printf("Original gap pos=%d len=%d\n", gap_iter.original_gap_position, gap_iter.original_gap_length);
            
            // Allocate space for index
            // Allocate counts[0][0][0] if not already allocated
            if (!counts[0][0][0]) {
                counts[0][0][0] = calloc(1, sizeof(*counts[0][0][0]));
         
            // Allocate counts[0][0][0][0] if not already allocated
            if (!counts[0][0][0][0]) {
                counts[0][0][0][0] = calloc(1, sizeof(*counts[0][0][0][0]));
            }
            // Allocate counts[0][0][0][0][0] if not already allocated
            if (!counts[0][0][0][0][0]) {
                counts[0][0][0][0][0] = calloc(1, sizeof(*counts[0][0][0][0][0]));
            }
            // Allocate counts[0][0][0][0][0][0] if not already allocated
            if (!counts[0][0][0][0][0][0]) {
                short int initial_indices = 100;
                counts[0][0][0][0][0][0] = calloc(initial_indices, sizeof(long int));
                counts[0][0][0][0][0][0][0] = initial_indices;  // size of allocation
                counts[0][0][0][0][0][0][1] = 1;                // number of indices found + 1 (as two slots used for allocation and last position)
            }
            }
            // Initialize the first index to indicate number of indices
            
            // Allocate for each gap pattern
            do {
                int len = gap_iter.added_length;
                int seq_space = 1 << (2 * len); // 4^len

                long int *index_array = counts[0][0][0][0][0][0];
                long int current_count = index_array[1];
                long int idx = pack_index(gap_iter.original_kmer_number, len, gap_iter.gap_pos[0], gap_iter.gap_pos[1], gap_iter.gap_pos[2]);

                if (index_array[1] >= index_array[0] - 2) {
                    long int number_of_indices = 2 * index_array[0];
                   // printf("\nReallocating to size %li", number_of_indices); fflush(stdout);
                   //printf("Array count before realloc: %ld\n", index_array[0]); fflush(stdout);
                    counts[0][0][0][0][0][0] = realloc(counts[0][0][0][0][0][0], number_of_indices * sizeof(long int));
                    index_array = counts[0][0][0][0][0][0];
                    index_array[0] = number_of_indices;
                    // printf("Array count after realloc: %ld, indices %ld\n", index_array[0], index_array[1]); fflush(stdout);
                }
                index_array[1]++;
                index_array[index_array[1]] = idx;

                //printf(" added packed index %ld at position %ld", idx, index_array[1]); fflush(stdout);
                
                // Dimension 4
                if (!counts[f][l][len]) {
                    counts[f][l][len] = calloc(len + 1, sizeof(*counts[f][l][len]));
                    if (!counts[f][l][len]) {
                        fprintf(stderr, "Failed to allocate counts[%d][%d][%d]\n", f, l, len);
                        exit(EXIT_FAILURE);
                    }
                    total_bytes += (len + 1) * sizeof(*counts[f][l][len]);
                }

                // Dimension 5
                if (!counts[f][l][len][gap_iter.gap_pos[0]]) {
                    counts[f][l][len][gap_iter.gap_pos[0]] = calloc(len + 1, sizeof(*counts[f][l][len][gap_iter.gap_pos[0]]));
                    if (!counts[f][l][len][gap_iter.gap_pos[0]]) {
                        fprintf(stderr, "Failed to allocate counts[%d][%d][%d][%d]\n", f, l, len, gap_iter.gap_pos[0]);
                        exit(EXIT_FAILURE);
                    }
                    total_bytes += (len + 1) * sizeof(*counts[f][l][len][gap_iter.gap_pos[0]]);
                }

                // Dimension 6
                if (!counts[f][l][len][gap_iter.gap_pos[0]][gap_iter.gap_pos[1]]) {
                    counts[f][l][len][gap_iter.gap_pos[0]][gap_iter.gap_pos[1]] =
                        calloc(len + 1, sizeof(*counts[f][l][len][gap_iter.gap_pos[0]][gap_iter.gap_pos[1]]));
                    if (!counts[f][l][len][gap_iter.gap_pos[0]][gap_iter.gap_pos[1]]) {
                        fprintf(stderr, "Failed to allocate counts[%d][%d][%d][%d][%d]\n",
                                f, l, len, gap_iter.gap_pos[0], gap_iter.gap_pos[1]);
                        exit(EXIT_FAILURE);
                    }
                    total_bytes += (len + 1) * sizeof(*counts[f][l][len][gap_iter.gap_pos[0]][gap_iter.gap_pos[1]]);
                }

                // Dimension 7
                if (!counts[f][l][len][gap_iter.gap_pos[0]][gap_iter.gap_pos[1]][gap_iter.gap_pos[2]]) {
                    counts[f][l][len][gap_iter.gap_pos[0]][gap_iter.gap_pos[1]][gap_iter.gap_pos[2]] =
                        calloc(seq_space, sizeof(*counts[f][l][len][gap_iter.gap_pos[0]][gap_iter.gap_pos[1]][gap_iter.gap_pos[2]]));
                    if (!counts[f][l][len][gap_iter.gap_pos[0]][gap_iter.gap_pos[1]][gap_iter.gap_pos[2]]) {
                        fprintf(stderr, "Failed to allocate counts[%d][%d][%d][%d][%d][%d]\n",
                                f, l, len, gap_iter.gap_pos[0], gap_iter.gap_pos[1], gap_iter.gap_pos[2]);
                        exit(EXIT_FAILURE);
                    }
                    total_bytes += seq_space * sizeof(*counts[f][l][len][gap_iter.gap_pos[0]][gap_iter.gap_pos[1]][gap_iter.gap_pos[2]]);
                }

                // printf("Allocating pattern: len=%d gaps=%d,%d,%d seq_space=%d\n",len,gap_iter.gap_pos[0],gap_iter.gap_pos[1],gap_iter.gap_pos[2],seq_space);
                
                // add original kmer counts to 7d matrix
                // if (f == 1 && !counts[f][l][0][0][0][0]) counts[f][l][0][0][0][0] = calloc (3, sizeof(long int));
                // if (f == 1) counts[f][l][0][0][0][0][0] = table->entries[l].signal_count;
                
            } while (move_gaps(&gap_iter));
        }
    }

    // printf("Total memory allocated: %ld bytes (%.2f MB)\n",total_bytes, total_bytes / (1024.0 * 1024.0));
    long int *index_array = counts[0][0][0][0][0][0];
    //printf("\nSorting"); fflush(stdout);
    return counts;
}



// Function to read local maxima table
LocalMaxTable* old_read_local_max_table(const char *filename) {
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        fprintf(stderr, "Could not open file %s\n", filename);
        return NULL;
    }
    
    LocalMaxTable *table = malloc(sizeof(LocalMaxTable));
    table->entries = NULL;
    table->num_entries = 0;
    
    char line[MAX_LINE_LEN];
    int capacity = 10;
    table->entries = malloc(capacity * sizeof(LocalMaxEntry));
    
    // Skip first two lines
    fgets(line, sizeof(line), fp);
    fgets(line, sizeof(line), fp);
    
    while (fgets(line, sizeof(line), fp)) {
        if (table->num_entries >= capacity) {
            capacity *= 2;
            table->entries = realloc(table->entries, capacity * sizeof(LocalMaxEntry));
        }
        
        char *token = strtok(line, " \t");
        int col = 0;
        
        while (token) {
            switch(col) {
                case 0:
                    table->entries[table->num_entries].sequence = strdup(token);
                    table->entries[table->num_entries].kmer_length = strlen(token);
                    break;
                case 2:
                    table->entries[table->num_entries].background_count = atoi(token);
                    break;
                case 3:
                    table->entries[table->num_entries].signal_count = atoi(token);
                    break;
            }
            token = strtok(NULL, " \t");
            col++;
        }
        table->num_entries++;
    }
    
    fclose(fp);
    return table;
}

// Function to check if character is a valid nucleotide (including IUPAC codes)
int is_nucleotide(char c) {
    // const char *valid_chars = "ACGTNRYMKWSBDHVacgtnrymkwsbdhv";
    const char *valid_chars = "ACGTNRYMKWSBDHVn";  // BASE, IUPAC OR LOWER CASE n. To prevent text from being misrecognized
    return strchr(valid_chars, c) != NULL;
}

// Function to check if string is a valid kmer
int is_valid_kmer(const char* sequence) {
    if (!sequence) return 0;
    
    for (const char *p = sequence; *p != '\0'; p++) {
        if (!is_nucleotide(*p)) return 0;
    }
    return 1;
}

// Function to read local maxima table
LocalMaxTable* read_local_max_table(const char *filename) {
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        fprintf(stderr, "Could not open file %s\n", filename);
        return NULL;
    }
    
    // Initialize table structure
    LocalMaxTable *table = malloc(sizeof(LocalMaxTable));
    if (!table) {
        fclose(fp);
        return NULL;
    }
    table->entries = NULL;
    table->num_entries = 0;
    
    // Initial allocation
    size_t capacity = 10;
    table->entries = malloc(capacity * sizeof(LocalMaxEntry));
    if (!table->entries) {
        free(table);
        fclose(fp);
        return NULL;
    }
    
    char line[MAX_LINE_LEN];
    while (fgets(line, sizeof(line), fp)) {
        char *first_token = strtok(line, " \t\n\r");
        if (!first_token || !is_valid_kmer(first_token)) {
            continue;
        }
        
        // Reallocate if needed
        if (table->num_entries >= capacity) {
            capacity *= 2;
            LocalMaxEntry *new_entries = realloc(table->entries, capacity * sizeof(LocalMaxEntry));
            if (!new_entries) {
                //free_local_max_table(table);
                fclose(fp);
                return NULL;
            }
            table->entries = new_entries;
        }
        
        // Process columns
        table->entries[table->num_entries].sequence = strdup(first_token);
        table->entries[table->num_entries].kmer_length = strlen(first_token);
        
        char *token;
        int col = 1;  // Start from column 1 since we already processed column 0
        
        while ((token = strtok(NULL, " \t\n\r"))) {
            if (col == 2) {
                table->entries[table->num_entries].background_count = atoi(token);
            } else if (col == 3) {
                table->entries[table->num_entries].signal_count = atoi(token);
            }
            col++;
        }
        
        table->num_entries++;
    }
    
    fclose(fp);
    return table;
}

// Helper function to convert sequence to index
long int sequence_to_index(const char* sequence) {
    long int index = 0;
    printf("Converting sequence '%s' to index:\n", sequence);
    for(int i = 0; sequence[i]; i++) {
        if(sequence[i] != 'n') {
            index = index << 2;
            switch(sequence[i]) {
                case 'A': printf("  A -> 00\n"); index += 0; break;
                case 'C': printf("  C -> 01\n"); index += 1; break;
                case 'G': printf("  G -> 10\n"); index += 2; break;
                case 'T': printf("  T -> 11\n"); index += 3; break;
            }
        } else {
            printf("  Skipping gap\n");
        }
    }
    printf("Final index: %ld\n", index);
    fflush(stdout);
    return index;
}

void find_extended_kmers(
   const char* sequence,
   int seq_len,
   const LocalMaxEntry* local_max,
   int local_max_idx,
   int file_num,
   long int******* counts,
   int max_length)
{
   int local_max_len = strlen(local_max->sequence);
   char *extracted_bases = malloc(MAX_KMER_LEN + 1);
   char *rc_sequence = malloc(seq_len + 1);

   // Build reverse complement
   for(int i = 0; i < seq_len; i++) {
       switch(sequence[seq_len - 1 - i]) {
           case 'A': rc_sequence[i] = 'T'; break;
           case 'T': rc_sequence[i] = 'A'; break;
           case 'G': rc_sequence[i] = 'C'; break;
           case 'C': rc_sequence[i] = 'G'; break;
           case 'n': rc_sequence[i] = 'n'; break;
           default : rc_sequence[i] = sequence[seq_len - 1 - i];
       }
   }
   rc_sequence[seq_len] = '\0';

   const char* sequences[2] = {sequence, rc_sequence};

   // printf("\nProcessing local max: %s (length %d)\n", local_max->sequence, local_max_len); fflush(stdout);

   // Find gap position and length in local max
   int gap_position = 0, gap_length = 0;
   for(int i = 0; i < local_max_len; i++) {
       if(local_max->sequence[i] == 'n') {
           if(gap_length == 0) gap_position = i;
           gap_length++;
       }
   }

   for(int seq_idx = 0; seq_idx < 2; seq_idx++) {
       const char* current_seq = sequences[seq_idx];

       for(int pos = 0; pos < seq_len - local_max_len + 1; pos++) {
           // Check for match with the local max sequence
           int match = 1;
           for(int i = 0; i < local_max_len; i++) {
               if(local_max->sequence[i] != 'n' &&
                  current_seq[pos + i] != local_max->sequence[i]) {
                   match = 0;
                   break;
               }
           }

           if(match) {
               // printf("\nFound match at position %d in %s sequence\n", pos, seq_idx == 0 ? "forward" : "reverse");

               // Initialize gap iterator for this match
               SingleAndTripleGapKmer gap_iter;
               gap_iter.original_kmer_length = local_max_len-gap_length;
               gap_iter.original_gap_position = gap_position;
               gap_iter.original_gap_length = gap_length;
               gap_iter.added_length = 0;
               gap_iter.max_length = max_length;
               for(int i = 0; i < 3; i++) {
                   gap_iter.gap_pos[i] = 0;
                   gap_iter.gap_len[i] = 0;
               }
               // runs iterator once to initialize and get rid of added length 0
               move_gaps(&gap_iter);
               
               // Try each possible gap pattern
               do {
                   int start_pos = pos - gap_iter.added_length;
                   int end_pos = pos + local_max_len + gap_iter.added_length;

                   {
                       //printf("\nTrying pattern: len=%d pos:len at %d:%d,%d:%d,%d:%d\n",gap_iter.added_length, gap_iter.gap_pos[0], gap_iter.gap_len[0], gap_iter.gap_pos[1], gap_iter.gap_len[1], gap_iter.gap_pos[2], gap_iter.gap_len[2]);
                       //printf("Match context in %ibp seq: ", seq_len);
                       // for(int i = start_pos; i < end_pos; i++) {
                       //   if(i == pos) printf("[");
                       //    if(i >= 0 && i < seq_len) printf("%c", current_seq[i]);
                       //    if(i == pos + local_max_len - 1) printf("]");
                       //}
                       // printf("\n");
                       //if (gap_iter.added_length == 2) exit(0);
                       
                       // Extract bases according to current gap pattern
                       int base_count = 0;
                       
                       int read_pos = pos - gap_iter.gap_pos[0];
                       
                       // continues if beginning of extracted kmer is outside bounds
                       if (read_pos < 0) continue;
                       // continues if end of extracted kmer is outside bounds
                       if (pos + gap_iter.original_kmer_length + gap_iter.original_gap_length + (gap_iter.added_length - gap_iter.gap_pos[2]) > seq_len) continue;
                   
                       for (int i = 0; i < gap_iter.gap_pos[0]; i++) {
                           // printf ("%c", current_seq[read_pos]);
                           extracted_bases[base_count++] = current_seq[read_pos++];
                       }
                       // printf ("-");
                   
                       read_pos += gap_iter.gap_len[0];
                       for (int i = gap_iter.gap_pos[0]; i < gap_iter.gap_pos[1]; i++) {
                           // printf ("%c", current_seq[read_pos]);
                           extracted_bases[base_count++] = current_seq[read_pos++];
                       }
                       // printf ("-");
                   
                       read_pos += gap_iter.gap_len[1];
                       for (int i = gap_iter.gap_pos[1]; i < gap_iter.gap_pos[2]; i++) {
                           // printf ("%c", current_seq[read_pos]);
                           extracted_bases[base_count++] = current_seq[read_pos++];
                       }
                       // printf ("-");
                       
                       read_pos += gap_iter.gap_len[2];
                       //if (gap_iter.gap_len[2] < 0) printf("\negative value for gap 3");
                                       for (; base_count < gap_iter.added_length; ) {
                                           //printf ("%c", current_seq[read_pos]);
                                           extracted_bases[base_count++] = current_seq[read_pos++];
                                           }
                       // printf ("\n");
                       
                       extracted_bases[base_count] = '\0';


                       // Convert to index
                       long int seq_index = 0;
                       for(int i = 0; i < base_count; i++) {
                           seq_index = seq_index << 2;
                           switch(extracted_bases[i]) {
                               case 'A': seq_index += 0; break;
                               case 'C': seq_index += 1; break;
                               case 'G': seq_index += 2; break;
                               case 'T': seq_index += 3; break;
                           }
                       }
                       // printf("\nSequence index: %ld\n", seq_index);

                       // if (gap_iter.added_length == 1) { printf("For kmer %s extracted: %s index %ld gaps: %d,%d,%d len: %d,%d,%d base_count: %d\n", local_max->sequence, extracted_bases, seq_index, gap_iter.gap_pos[0], gap_iter.gap_pos[1], gap_iter.gap_pos[2], gap_iter.gap_len[0], gap_iter.gap_len[1], gap_iter.gap_len[2], base_count);}
                       
                       //printf("Incrementing: len=%d gaps=%d,%d,%d\n", gap_iter.added_length, gap_iter.gap_pos[0], gap_iter.gap_pos[1], gap_iter.gap_pos[2]);
                       // fflush(stdout);
                       
                       if (!counts[file_num][local_max_idx][gap_iter.added_length]
                           [gap_iter.gap_pos[0]][gap_iter.gap_pos[1]]
                           [gap_iter.gap_pos[2]]) printf ("\nMemory not allocated!");
                       
                       else counts[file_num][local_max_idx][gap_iter.added_length]
                             [gap_iter.gap_pos[0]][gap_iter.gap_pos[1]]
                             [gap_iter.gap_pos[2]][seq_index]++;
                       
                       // printf("\nNext step"); fflush(stdout);
                       
                   }
               } while(move_gaps(&gap_iter));
           }
       }
   }

   free(extracted_bases);
   free(rc_sequence);
}

void process_sequence_file(const char* filename,
                         LocalMaxTable* table,
                         long int******* counts,
                         int file_num,
                         int max_extension)
{
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        fprintf(stderr, "Could not open file %s\n", filename);
        return;
    }
    
    char *sequence = malloc(MAX_LINE_LEN);
    while (fgets(sequence, MAX_LINE_LEN, fp)) {
        sequence[strcspn(sequence, "\n")] = 0;
        int seq_len = strlen(sequence);
        
        for (int i = 0; i < table->num_entries; i++) {
            find_extended_kmers(sequence, seq_len, &table->entries[i], i,
                              file_num, counts, max_extension);
        }
    }
    
    free(sequence);
    fclose(fp);
}

// Calculate information content (2-entropy) from array of 4 frequencies
double calculate_position_ic(double freqs[4]) {
   double ic = 0.0;
   for(int i = 0; i < 4; i++) {
       if(freqs[i] > 0) {
           ic -= freqs[i] * log2(freqs[i]);
       }
   }
   return 2.0 - ic;
}

double calculate_position_variance(const uint64_t counts[4])
{
    // Number of data points
    const int n = 4;

    // 1) Compute the mean of all four values
    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        sum += (double)counts[i];
    }
    double mean = sum / (double)n;

    // 2) Sum up the squared differences from the mean
    double var_sum = 0.0;
    for (int i = 0; i < n; i++) {
        double diff = (double)counts[i] - mean;
        var_sum += diff * diff;
    }

    // 3) Sample variance = var_sum / (n - 1)
    //    Sample standard deviation = sqrt(sample variance)
    if (n > 1) {
        double sample_variance = var_sum / (double)(n);
        return sqrt(sample_variance);
    } else {
        // In a degenerate case (n <= 1), just return 0
        return 0.0;
    }
}

short int print_kmer_statistics(long int *count_array, int file_number, int kmer_length, long int kmer) {
   // Original kmer count
    //return;
    
   uint64_t orig_count = count_array[kmer];
      
   // 2. Calculate IC and position-specific variance
   double total_ic = 0.0;
   double total_var = 0.0;
   int valid_positions = 0;
   int hamming_max = 1;
    long int kmer1_count = count_array[kmer];
    
   // For each position
   for(int pos = 0; pos < kmer_length; pos++) {
       // Get counts for consensus and mutations at this position
       uint64_t pos_counts[4] = {0};  // A,C,G,T counts
       
       int shift = 2 * pos;
       
       // Get counts for consensus and all mutations at this position
       uint64_t base_mask = ~(3ULL << shift);
       for(int base = 0; base < 4; base++) {
           uint64_t mut_kmer = (kmer & base_mask) | ((uint64_t)base << shift);
           pos_counts[base] = count_array[mut_kmer];
           if (pos_counts[base] > kmer1_count) hamming_max = 0;
       }
       
       // Calculate frequencies for IC
       double total = 0.0;
       double freqs[4] = {0.0};
       for(int i = 0; i < 4; i++) {
           total += (double)pos_counts[i];
       }
       if(total > 0.0) {
           for(int i = 0; i < 4; i++) {
               freqs[i] = (double)pos_counts[i]/total;
           }
       }
       
       total_ic += calculate_position_ic(freqs);
       total_var += calculate_position_variance(pos_counts);
       valid_positions++;
   }
   
   double avg_variance = valid_positions > 0 ? total_var/(double)valid_positions : 0.0;
   double exp_variance = sqrt((double)orig_count);  // Expected Poisson variance
   if (exp_variance == 0) exp_variance = 0.001;
    
   // Print tab-separated statistics
printf("\t%.3lf\t%.2lf%%\t%.2lf%%\t%.2lf%%",
               total_ic,
               100 * avg_variance/orig_count,
               100 * exp_variance/orig_count,
               100 * avg_variance/orig_count - 100 * exp_variance/orig_count);
 
    return hamming_max;
}

int main(int argc, char *argv[]) {
    if (argc < 4 || argc > 8) {
        fprintf(stderr, "Usage: %s background_seq signal_seq local_max_table [max_extension] [min count for print] [length_diff_cutoff] [r or robust if robust Localmax is required]\n", argv[0]);
        return 1;
    }
    
    const char *background_file = argv[1];
    const char *signal_file = argv[2];
    const char *local_max_file = argv[3];
    int max_extension = (argc >= 5) ? atoi(argv[4]) : 2;
    int print_count_cutoff = (argc >= 6) ? atoi(argv[5]) : 100;
    double length_diff_cutoff = (argc >= 7) ? atof(argv[6]) : 0.35;
    short int robust = (argc >= 8) ? argv[7][0] == 'r' : 0;
    
    //printf ("\nrobust = %i", robust); fflush(stdout);
    if(max_extension > MAX_ADDED_LEN) {
        printf("\nToo long extension, max set to %i \nEdit source code parameter MAX_ADDED_LEN or extend less\n", MAX_ADDED_LEN);
        exit(1);
    }
    
    short int scale_factor = 1;
    short int background_subtraction = 0;
    
    LocalMaxTable *table = read_local_max_table(local_max_file);
    if (!table) return 1;
    
    GenerateMask (); // generates mask for bit operations
    
    long int******* counts = allocate_counts(2, table, max_extension); // allocates memory for counts
    
    // processes sequence files to count all triplegap kmers
    process_sequence_file(background_file, table, counts, 0, max_extension);
    process_sequence_file(signal_file, table, counts, 1, max_extension);
    
    // Print results for each local max
    // printf("\nExtended k-mers found:\n");
    printf("\nFull_kmer\tBackground count\tSignal count\tOriginal_count\tOriginal_localmax\tExtender");
    
    for (int l = 0; l < table->num_entries; l++) {
        
        MergedKmer merged_kmer;
        // Initialize gap iterator
        SingleAndTripleGapKmer gap_iter;
        Set_original_gap_info_and_init(table->entries[l].sequence, &gap_iter, max_extension, l);
        
        /* Print original kmer if no kmer that is longer by 1 is localmax
        gap_iter.triplegap_kmer_value = 0; // initialize value
        merged_kmer.kmer_value = gap_iter.original_kmer_value;
        merged_kmer.kmer_length = gap_iter.original_kmer_length;
        merged_kmer.gap_position = gap_iter.original_gap_position;
        merged_kmer.gap_length = gap_iter.original_gap_length;
        SingleAndTripleGapKmer *temp_triplegap_kmer = malloc(sizeof(SingleAndTripleGapKmer));
        memcpy(temp_triplegap_kmer, &gap_iter, sizeof(SingleAndTripleGapKmer));
        //temp_triplegap_kmer->debug = 1;
        //if (Mergedkmer_Localmax(counts, temp_triplegap_kmer, &merged_kmer, background_subtraction, scale_factor, robust, length_diff_cutoff)) printf("\n%s\t%i\t%i\t%i\t%s\tNo_change\tNA\tNA\tNA\tNA\tHamming_max\tHigh\tLocalmax", table->entries[l].sequence, table->entries[l].background_count, table->entries[l].signal_count, table->entries[l].signal_count, table->entries[l].sequence);
        //temp_triplegap_kmer->debug = 0;
        free(temp_triplegap_kmer); */
        
        
        // sort allocation array so that localmax in print can find which kmers are allocated (counted)
        qsort(counts[0][0][0][0][0][0] + 2, counts[0][0][0][0][0][0][1] - 1, sizeof(long int), compare_longs);
        //printf("\nAllocation complete"); fflush(stdout);
        
        short int by1_extended_kmer_has_higher_counts = 0;
        
        // For each valid gap pattern
        do {
            long int total_number_of_kmers = 1 << (2 * gap_iter.added_length);
            for (long int kmer = 0; kmer < total_number_of_kmers; kmer++)
            {
            long int background_count = counts[0][l][gap_iter.added_length]
                                     [gap_iter.gap_pos[0]]
                                     [gap_iter.gap_pos[1]]
                                     [gap_iter.gap_pos[2]][kmer];
            long int signal_count = counts[1][l][gap_iter.added_length]
                                 [gap_iter.gap_pos[0]]
                                 [gap_iter.gap_pos[1]]
                                 [gap_iter.gap_pos[2]][kmer];
            
                if (signal_count > print_count_cutoff) {
                    gap_iter.triplegap_kmer_value = kmer;
                    printf ("\n");
                    Merge_kmers(&gap_iter, &merged_kmer, 1);
                    //printf(" gap pos:len %i:%i,%i:%i,%i:%i ", gap_iter.gap_pos[0], gap_iter.gap_len[0],gap_iter.gap_pos[1], gap_iter.gap_len[1], gap_iter.gap_pos[2], gap_iter.gap_len[2]);
                    //Align_to_identify_triplegapkmer(&gap_iter, &merged_kmer);
                    printf ("\t");
                    // Test_Compute_triplegap_kmer(&gap_iter, &merged_kmer);
                    printf("\t%ld\t%li\t%i\t", background_count, signal_count, table->entries[l].signal_count);
                    Kmerprint(gap_iter.original_kmer_value, gap_iter.original_kmer_length, gap_iter.original_gap_position, gap_iter.original_gap_length);
                    printf ("\t");
                    Triplegap_kmerprint(&gap_iter);
                    if (print_kmer_statistics(counts[1][l][gap_iter.added_length]
                                          [gap_iter.gap_pos[0]]
                                          [gap_iter.gap_pos[1]]
                                              [gap_iter.gap_pos[2]], 1, gap_iter.added_length, gap_iter.triplegap_kmer_value)) printf ("\tHamming_max");
                    else printf("\t");
                    short int higher_than_original = signal_count > (table->entries[l].signal_count * pow(length_diff_cutoff, gap_iter.added_length));
                    if (higher_than_original) {
                        printf("\tHigh");
                        if (gap_iter.added_length == 1) {
                            by1_extended_kmer_has_higher_counts = 1;
                            //printf("\tHigher_than_original");
                        } // at least one that has higher counts
                    }
                    SingleAndTripleGapKmer *temp_triplegap_kmer = malloc(sizeof(SingleAndTripleGapKmer));
                    memcpy(temp_triplegap_kmer, &gap_iter, sizeof(SingleAndTripleGapKmer));
                    if (higher_than_original && Mergedkmer_Localmax(counts, temp_triplegap_kmer, &merged_kmer, background_subtraction, scale_factor, robust, length_diff_cutoff)) printf("\tLocalmax");
                    free(temp_triplegap_kmer);

                }
                //printf ("\n");
            }
            
        } while (move_gaps(&gap_iter));
        
        // original kmer becomes Localmax if no kmer extended by 1 has higher count
        if (by1_extended_kmer_has_higher_counts == 0) printf("\n%s\t%i\t%i\t%i\t%s\tNo_change\tNA\tNA\tNA\tNA\tHamming_max\tHigh\tLocalmax", table->entries[l].sequence, table->entries[l].background_count, table->entries[l].signal_count, table->entries[l].signal_count, table->entries[l].sequence);
    }
    
    // Cleanup
    for (int i = 0; i < table->num_entries; i++) {
        free(table->entries[i].sequence);
    }
    free(table->entries);
    free(table);
    
    printf ("\n");
    
    return 0;
}
