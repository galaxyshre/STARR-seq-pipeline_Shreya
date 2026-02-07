#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <pthread.h>
#include <inttypes.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <stdatomic.h>

#define MAX(x, y) ((x) > (y) ? (x) : (y))

/* VERSION 0.2, 13 Jan 2025   part of AUTOSEED code modified to enable bitstream input (sequences of arbitrary length)   */
/* bitstream and multithreaded code written with assistance from Claude 3.5 Sonnet and ChatGPT4 and o1                   */
/* Debugged manually by J Taipale as bit shifting is not their forte                                                     */
/* Compile with clang -march=native -ffast-math -O3 -o localmax_motif localmax_motif.c                             */

/* GLOBAL VARIABLES */
uint64_t mask_ULL[42][42];   /* PRIMARY mask_ULL FOR EACH SEPARATE NUCLEOTIDE */
uint64_t lowmask_ULL[42];    /* LOW MASK FOR EACH KMER */
uint64_t highmask_ULL[42];   /* HIGH MASK FOR EACH KMER */
short int Nlength = 32;              /* LENGTH OF MASK */

#define NUM_FILES 2
#define MAX_GAP_LENGTH 10
#define KMER_VARIATION 1
#define BUFFER_SIZE 1024
#define MAX_KMER_LEN 32  // Adjust if needed for your use case
#define MAX_SEQ_LEN 500  // Adjust if needed for your use case
#define MAX_WIDTH_OF_PWM 40
#define FLANK_WIDTH 5
#define SHIFTED_GAP_POSITIONS 1    // count also gap positions that are shifted from center

// STRUCTURES

typedef _Atomic size_t atomic_size_t;

// LineBuffer structure to include 5D results
typedef struct {
    uint64_t***** results;  // 5D array [file][kmer_len][gap_pos][gap_len][kmer]
    uint64_t* bitstream;
    uint64_t* temp_bitstream;
    int max_seq_length;
    int file_index;
    char finished;
    pthread_mutex_t mutex;
    pthread_cond_t not_full;
    pthread_cond_t not_empty;
    char** sequences;
    size_t* seq_lengths;
    int write_pos;
    int read_pos;
    int count;
    int shortest_kmer;     // Add these
    int longest_kmer;      // two fields
    FILE* file;
} LineBuffer;


typedef struct PWMConsumerArgs {
    LineBuffer* buffer;
    struct normalized_pwm* signal_pwm;        // Input PWM for finding matches
    struct normalized_pwm* background_pwm;     // Input background PWM for corrections
    struct count_pwm* signal_output_pwm;      // Output PWM for signal counts
    struct count_pwm* background_output_pwm;   // Output PWM for background counts
    double cutoff;
    _Atomic(size_t) *total_reads;    // Use _Atomic qualifier
    _Atomic(size_t) *matched_reads;
    int num_files;
    // K-mer counting additions
    int kmer_length;                 // Length of k-mers to count
    double lambda;
    uint64_t** kmer_counts;          // Array to store k-mer counts for all files
} PWMConsumerArgs;


struct match {short int *position; double *score;};
short int match_init (struct match *i, short int width)
{
short int maximum_width = Nlength+10;
short int counter;
(*i).position = malloc(sizeof(short int) * maximum_width + 5);
(*i).score = malloc(sizeof(double) * maximum_width + 5);
for (counter = 0; counter < maximum_width; counter++)
{
(*i).position[counter] = 0;
(*i).score[counter] = 0;
}
return(0);
}

/* COUNT PWM */
struct count_pwm {char *name; short int width; long int max_counts; double **incidence; double enrichment;};
short int count_pwm_clear (struct count_pwm *i, char *name, short int width, double initial_value)
{
short int maximum_width = MAX_WIDTH_OF_PWM;
short int counter;
short int counter2;
strcpy ((*i).name, name);
(*i).enrichment = 0;
(*i).width = width;
(*i).max_counts = initial_value;
for (counter = 0; counter < 5; counter++)
{
for (counter2 = 0; counter2 < maximum_width; counter2++) (*i).incidence[counter][counter2] = initial_value;
}
return(0);
}
short int count_pwm_init (struct count_pwm *i, char *name, short int width, double initial_value)
{
short int maximum_width = MAX_WIDTH_OF_PWM;
short int counter;
short int counter2;
(*i).name = malloc(1000);
strcpy ((*i).name, name);
(*i).width = width;
(*i).max_counts = initial_value;
(*i).incidence = malloc(sizeof(double *) * 5 + 5);
for (counter = 0; counter < 5; counter++)
{
(*i).incidence[counter] = malloc(sizeof(double) * maximum_width + 5);
for (counter2 = 0; counter2 < maximum_width; counter2++) (*i).incidence[counter][counter2] = initial_value;
}
return(0);
}
short int count_pwm_free (struct count_pwm *i)
{
    short int counter;
    free((*i).name);
    for (counter = 0; counter < 5; counter++) free((*i).incidence[counter]);
    free((*i).incidence);
    return(0);
}

/* NORMALIZED PWM */
struct normalized_pwm {char *name; char *seed; short int width; long int max_counts; double *information_content; short int *original_position; double *position_score; long int *total_counts_for_column; double **fraction; double enrichment; short int negative_values_allowed;};
short int normalized_pwm_init (struct normalized_pwm *i, char *name, short int width, double initial_value)
{
short int maximum_width = MAX_WIDTH_OF_PWM;
short int counter;
short int counter2;
(*i).enrichment = 0;
(*i).negative_values_allowed = 0;
(*i).name = malloc(100);
strcpy ((*i).name, name);
(*i).seed = malloc(1000);
strcpy ((*i).seed, "UNKNOWN");
(*i).width = width;
(*i).max_counts = initial_value;
(*i).fraction = malloc(sizeof(double *) * 5 + 5);
(*i).information_content = malloc(sizeof(double) * maximum_width + 5);
(*i).position_score = malloc(sizeof(double) * maximum_width + 5);
(*i).original_position = malloc(sizeof(short int) * maximum_width + 5);
(*i).total_counts_for_column = malloc(sizeof(long int) * maximum_width + 5);

for (counter = 0; counter < 5; counter++)
{
(*i).fraction[counter] = malloc(sizeof(double) * maximum_width + 5);
for (counter2 = 0; counter2 < maximum_width; counter2++) (*i).fraction[counter][counter2] = initial_value;
}
for (counter2 = 0; counter2 < maximum_width; counter2++)
{
(*i).information_content[counter2] = 0;
(*i).position_score[counter2] = 0;
(*i).original_position[counter2] = counter2;
(*i).total_counts_for_column[counter2] = 0;
}
return(0);
}
short int normalized_pwm_free (struct normalized_pwm *i)
{
short int counter;
free((*i).name);
free((*i).information_content);
free((*i).position_score);
free((*i).total_counts_for_column);
for (counter = 0; counter < 5; counter++) free((*i).fraction[counter]);
free((*i).fraction);
return(0);
}

/* SUBROUTINE THAT GENERATES PWM FROM IUPAC */
short int Iupac_to_pwm(struct normalized_pwm *n, char *searchstring)
{
short int counter;
short int pwm_position;
short int current_match_position;
short int nucleotide_value;
(*n).width =  strlen(searchstring);

char **canbe;
canbe = malloc(sizeof(char *) * 4 + 5);
for (counter = 0; counter < 4; counter++) canbe[counter] = malloc(200);
strcpy (canbe[0], "100010101001111");
strcpy (canbe[1], "010001100110111");
strcpy (canbe[2], "001010010111011");
strcpy (canbe[3], "000101011011101");

char *conversion_string;
conversion_string = malloc(sizeof(char)*3 + 5);
conversion_string[0] = '\0';
conversion_string[1] = '\0';

char *nucleotide_iupac = "ACGTRYMKWSBDHVN";
short int iupac_length = 15;
    
/* BUILDS PWM */
for(pwm_position = 0; pwm_position < (*n).width; pwm_position++)
{
// substitutes to allow lower case n
    if (searchstring[pwm_position] == 'n') searchstring[pwm_position] = 'N';
    
//printf("\n");
for(current_match_position = 0; (current_match_position < iupac_length) & (searchstring[pwm_position] != nucleotide_iupac[current_match_position]); current_match_position++)
;
    //if(current_match_position == iupac_length) {strcat(seed_story, "\n** SEED ERROR: defective IUPAC"); printf("%s", seed_story); fflush(stdout); exit(1);}
for (nucleotide_value = 0; nucleotide_value < 4; nucleotide_value++)
{
conversion_string[0] = canbe[nucleotide_value][current_match_position];
(*n).fraction[nucleotide_value][pwm_position] = atof(conversion_string) - 1;
}
}

for (counter = 0; counter < 4; counter++) free(canbe[counter]);
free(conversion_string);
return(0);
}

// Print count motif
short int PrintMotif(struct count_pwm *multinomial_motif, short int multinomial, int seed_length) {
    //printf("\nDEBUG: count_pwm structure: width=%d\n", multinomial_motif->width);
    // printf("\nMultinomial %i motif counts, attempting to print all positions:\n", multinomial);
    // printf("Debug: PWM address in print: %p\n", (void*)multinomial_motif);
    // Print positions header
    printf("\n");
    for (int i = 1; i < multinomial_motif->width; i++) {
        printf("\t%d", i);
    }
    printf("\n");
    
    // Print each base as a row
    const char bases[] = {'A', 'C', 'G', 'T'};
    for (int base = 0; base < 4; base++) {
        //printf("%c  ", bases[base]);
        for (int i = 1; i < multinomial_motif->width; i++) {
            printf("\t%.0f", multinomial_motif->incidence[base][i]);
        }
        printf("\n");
    }
    return 0;
}

// Print count motif
short int PrintNormalizedMotif(struct normalized_pwm *motif) {
    // Print positions header
    printf("\n   ");
    for (int i = 0; i < motif->width; i++) {
        printf("\t%d", i);
    }
    printf("\n");
    
    // Print each base as a row
    const char bases[] = {'A', 'C', 'G', 'T'};
    for (int base = 0; base < 4; base++) {
        printf("%c  ", bases[base]);
        for (int i = 0; i < motif->width; i++) {
            printf("\t%.0f", motif->fraction[base][i]);
        }
        printf("\n");
    }
    return 0;
}

short int Findpwmmatch_from_bitstream(struct normalized_pwm *p, double cut_off,
                                    uint64_t *bitstream, size_t sequence_length,
                                    struct match *match) {
   
    //PrintNormalizedMotif(p);
    //printf("\nCutoff: %f", cut_off);
    //printf("\nSequence length: %li", sequence_length);
    //printf("\nSeed motif width: %hi", (*p).width);
    
    uint64_t main_reg = 0;
    uint64_t source_reg = 0;
    const int shift_amount = 2;
    short int number_of_matches = 0;
    short int nucleotide;
    size_t current_position = 0;
    double score;
    signed short int pwm_position;
    size_t last_source_block_moved = 0;
    size_t current_block = 0;
    
    size_t bitstream_length = (sequence_length * 2 + 63) / 64;
    uint64_t *current_bitstream = bitstream;

    // Load initial blocks
    if (bitstream_length > 0) {
        main_reg = *current_bitstream++;
        bitstream_length--;
    }
    if (bitstream_length > 0) {
        source_reg = *current_bitstream++;
        bitstream_length--;
    }

    size_t total_blocks = (sequence_length + 31) / 32;
    size_t total_number_of_blocks_left_to_move_to_source = total_blocks - 2;
    size_t remaining_bits_in_main = 64;
    
    // Process all possible positions
    while (current_position + (*p).width <= sequence_length) {
        // Calculate score for current position
        for (score = 0, pwm_position = (*p).width-1; pwm_position >= 0; pwm_position--) {
            nucleotide = (main_reg & mask_ULL[1][pwm_position]) >> (2 * pwm_position);
            score += (*p).fraction[nucleotide][(*p).width-1-pwm_position];
        }
        //printf("\nScore %f at position: %lu", score, sequence_length - current_position - (*p).width);
        // Record match if score exceeds cutoff
        if (score >= cut_off) {
            number_of_matches++;
            (*match).score[number_of_matches] = score;
            (*match).position[number_of_matches] = sequence_length - current_position - (*p).width;
            //printf("\nMatch!! Score %f at position: %hi",score, (*match).position[number_of_matches]);
        }
        
        // Shift bits
        uint64_t new_bits = (source_reg & 0x3) << (64 - shift_amount);
        main_reg = (main_reg >> shift_amount) | new_bits;
        source_reg >>= shift_amount;
        remaining_bits_in_main -= 2;
        current_position++;

        // Handle last block transition
        if (last_source_block_moved == 0 && remaining_bits_in_main == 0 &&
            total_number_of_blocks_left_to_move_to_source == 0) {
            remaining_bits_in_main = (sequence_length % 32) * 2;
            last_source_block_moved = 1;
        }
        
        // Load next block if needed
        if (last_source_block_moved == 0 && remaining_bits_in_main == 0 &&
            total_number_of_blocks_left_to_move_to_source > 0) {
            source_reg = *current_bitstream++;
            total_number_of_blocks_left_to_move_to_source--;
            remaining_bits_in_main = 64;
            current_block++;
        }
    }
    
    (*match).position[0] = number_of_matches;
    return number_of_matches;
}

short int Multinomial_add_to_pwm_from_bitstream(struct count_pwm *p, struct normalized_pwm *qp,
                                              short int match_position, double score, double cut_off,
                                              uint64_t *bitstream, uint64_t *temp_bitstream, size_t sequence_length,
                                              struct normalized_pwm *background_pwm) {
    
    //printf("\nDebugging add_to_pwm: seq_length=%zu, match_position=%d, pwm width=%d\n", sequence_length, match_position, (*p).width);
    
    short int pwm_position;
    short int seed_position;
    short int background_position;
    short int query_position = 0;
    short int in_seed;
    uint64_t main_reg = 0;
    uint64_t source_reg = 0;
    const int shift_amount = 2;
    uint64_t nucleotide;
    size_t last_source_block_moved = 0;
    size_t current_block = 0;

    __atomic_add_fetch(&((*p).max_counts), 1, __ATOMIC_SEQ_CST);
    //(*p).max_counts++;
    
    size_t bitstream_length = (sequence_length * 2 + 63) / 64;
    for(int i = 0; i < bitstream_length; i++) temp_bitstream[i] = bitstream[i]; // copies bitstream
    uint64_t *current_bitstream = temp_bitstream;
    
    // Load initial blocks
    if (bitstream_length > 0) {
        main_reg = *current_bitstream++;
        bitstream_length--;
    }
    if (bitstream_length > 0) {
        source_reg = *current_bitstream++;
        bitstream_length--;
    }

    size_t total_blocks = (sequence_length + 31) / 32;  // How many blocks we should have
    size_t total_number_of_blocks_left_to_move_to_source = total_blocks - 2;
    size_t remaining_bits_in_main = 64;
    
    signed short int shifted_match_position = sequence_length - match_position - (*qp).width - 1;
    short int offset = 0;
    if (sequence_length - match_position <= FLANK_WIDTH + (*qp).width) offset = (*qp).width + FLANK_WIDTH - (sequence_length - match_position);
    
    pwm_position = 2 * FLANK_WIDTH + (*qp).width - offset;
    //printf("\nOffset: %i", offset);
    
    for (; shifted_match_position + (*qp).width + FLANK_WIDTH >= 0 && sequence_length > 1; shifted_match_position--, sequence_length--)
    {

        // shift collect bases from -10 until +10 bp flanking the seed match
        if (shifted_match_position < FLANK_WIDTH)
        {
            nucleotide = main_reg & 0x3;
            seed_position = pwm_position - FLANK_WIDTH;
            if (seed_position > 0 && seed_position <= (*qp).width) in_seed = 1;
            else {in_seed = 0; seed_position = 0;} // to prevent segfault
            
            // EXCLUDES CONSENSUS IF THERE IS MISMATCH
            if (in_seed == 1 && score <= cut_off + 1 && (*qp).fraction[nucleotide][seed_position-1] == 0)
            {
                // printf("\nExcluding nucleotide %c from PWM position %i, shifted match %i, seed position %i", "ACGT"[nucleotide], pwm_position, shifted_match_position, seed_position);
                // (*p).incidence[nucleotide][pwm_position]++;
                pwm_position--;
            }
            else
            {
                // printf("\nAdding nucleotide %c to PWM position %i, shifted match %i", "ACGT"[nucleotide], pwm_position, shifted_match_position);
                __atomic_add_fetch(&((*p).incidence[nucleotide][pwm_position]), 1, __ATOMIC_SEQ_CST);
                //(*p).incidence[nucleotide][pwm_position]++;
                pwm_position--;
            }
            }
        
            // Shift bits
            uint64_t new_bits = (source_reg & 0x3) << (64 - shift_amount);
            main_reg = (main_reg >> shift_amount) | new_bits;
            source_reg >>= shift_amount;
            remaining_bits_in_main -= 2;

            // Handle last block transition
            if (last_source_block_moved == 0 && remaining_bits_in_main == 0 &&
                total_number_of_blocks_left_to_move_to_source == 0) {
                remaining_bits_in_main = (sequence_length % 32) * 2;
                last_source_block_moved = 1;
            }
            
            // Load next block if needed
            if (last_source_block_moved == 0 && remaining_bits_in_main == 0 &&
                total_number_of_blocks_left_to_move_to_source > 0) {
                source_reg = *current_bitstream++;
                total_number_of_blocks_left_to_move_to_source--;
                remaining_bits_in_main = 64;
                current_block++;
            }

        
        }
    
    return(0);
}

void print_match_sequences(struct match *matches, uint64_t *bitstream, size_t sequence_length) {
    const size_t PRINT_LENGTH = 40;
    char sequence[PRINT_LENGTH + 1];  // +1 for null terminator
    sequence[PRINT_LENGTH] = '\0';
    
    uint64_t main_reg = 0;
    uint64_t source_reg = 0;
    const int shift_amount = 2;
    size_t last_source_block_moved = 0;
    
    // For each match
    for (int i = 1; i <= matches->position[0]; i++) {
        size_t match_pos = matches->position[i];
        size_t start_pos = 0;
        
        // Fill with 'x' first
        memset(sequence, 'x', PRINT_LENGTH);
        
        // Calculate actual sequence start and length
        if (match_pos > PRINT_LENGTH/2) {
            start_pos = match_pos - PRINT_LENGTH/2;
        }
        
        // Initialize bitstream reading
        size_t bitstream_length = (sequence_length * 2 + 63) / 64;
        uint64_t *current_bitstream = bitstream;
        
        // Load initial blocks
        if (bitstream_length > 0) {
            main_reg = *current_bitstream++;
            bitstream_length--;
        }
        if (bitstream_length > 0) {
            source_reg = *current_bitstream++;
            bitstream_length--;
        }

        size_t total_blocks = (sequence_length + 31) / 32;
        size_t total_number_of_blocks_left_to_move_to_source = total_blocks - 2;
        size_t remaining_bits_in_main = 64;
        
        // Skip to start position
        for (size_t pos = 0; pos < start_pos; pos++) {
            uint64_t new_bits = (source_reg & 0x3) << (64 - shift_amount);
            main_reg = (main_reg >> shift_amount) | new_bits;
            source_reg >>= shift_amount;
            remaining_bits_in_main -= 2;

            if (last_source_block_moved == 0 && remaining_bits_in_main == 0 &&
                total_number_of_blocks_left_to_move_to_source == 0) {
                remaining_bits_in_main = (sequence_length % 32) * 2;
                last_source_block_moved = 1;
            }
            
            if (last_source_block_moved == 0 && remaining_bits_in_main == 0 &&
                total_number_of_blocks_left_to_move_to_source > 0) {
                source_reg = *current_bitstream++;
                total_number_of_blocks_left_to_move_to_source--;
                remaining_bits_in_main = 64;
            }
        }
        
        // Extract sequence
        for (size_t pos = start_pos; pos < start_pos + PRINT_LENGTH && pos < sequence_length; pos++) {
            uint64_t nucleotide = main_reg & 0x3;
            sequence[pos - start_pos] = "ACGT"[nucleotide];
            
            // Shift to next position if not at end
            if (pos + 1 < start_pos + PRINT_LENGTH && pos + 1 < sequence_length) {
                uint64_t new_bits = (source_reg & 0x3) << (64 - shift_amount);
                main_reg = (main_reg >> shift_amount) | new_bits;
                source_reg >>= shift_amount;
                remaining_bits_in_main -= 2;

                if (last_source_block_moved == 0 && remaining_bits_in_main == 0 &&
                    total_number_of_blocks_left_to_move_to_source == 0) {
                    remaining_bits_in_main = (sequence_length % 32) * 2;
                    last_source_block_moved = 1;
                }
                
                if (last_source_block_moved == 0 && remaining_bits_in_main == 0 &&
                    total_number_of_blocks_left_to_move_to_source > 0) {
                    source_reg = *current_bitstream++;
                    total_number_of_blocks_left_to_move_to_source--;
                    remaining_bits_in_main = 64;
                }
            }
        }
        
        printf("%s\n", sequence);
    }
}


/* SUBROUTINE THAT PRINTS KMER (USES NO MASK) */
void Kmerprint (__uint128_t print_sequence_value, short int kmer_length, short int gap_position, short int gap_length)
{
    for (kmer_length--; kmer_length >= 0; kmer_length--, gap_position--)
    {   if (gap_position == 0) for(int i = 0; i < gap_length; i++) printf("n");
        printf("%c","ACGT"[(print_sequence_value >> (kmer_length * 2)) & 3]);
    }
}

// Comparison function for qsort
static int compare_counts(const void *a, const void *b) {
    return (*(long int*)b - *(long int*)a); // Sort in descending order
}

// smooths counts if smooth parameter is set to 2
long int smooth_kmer_counts(long int *****results,
                          short int file_number,
                          long int kmer,
                          short int kmer_length,
                          short int gap_position,
                          short int gap_length,
                          short int smooth) {
    
    // Return unsmoothed count if smooth is 0
    if (smooth != 2) {
        return results[file_number][kmer_length][gap_position][gap_length][kmer];
    }
    
    // Get original kmer count
    long int original_count = results[file_number][kmer_length][gap_position][gap_length][kmer];
    
    // Store neighbor counts (excluding original kmer)
    long int neighbor_counts[3 * kmer_length]; // Max neighbors: 3 options per position
    int count = 0;
    
    // Generate all Hamming distance 1 neighbors
    long int lowbit = 1;  // For A-C or G-T transitions
    long int highbit = 2; // For A-G or C-T transitions
    long int compared_kmer;
    
    for(short int position = 0; position < kmer_length; position++, lowbit <<= 2, highbit <<= 2) {
        // First substitution (A↔C, G↔T)
        compared_kmer = lowbit ^ kmer;
        neighbor_counts[count++] = results[file_number][kmer_length][gap_position][gap_length][compared_kmer];
        
        // Second substitution (A↔G, C↔T)
        compared_kmer = highbit ^ kmer;
        neighbor_counts[count++] = results[file_number][kmer_length][gap_position][gap_length][compared_kmer];
        
        // Third substitution (A↔T, C↔G)
        compared_kmer = lowbit ^ compared_kmer;
        neighbor_counts[count++] = results[file_number][kmer_length][gap_position][gap_length][compared_kmer];
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
    long int *****results,
    short int file_number,
    short int current_kmer_length,
    short int current_gap_position,
    short int current_gap_length,
    long int compared_kmer,
    long int kmer1_incidence,
    int background_subtraction,
    double scale_factor,
    double kmer_length_difference_cutoff,
    short int original_kmer_length)  // Added parameter
{
    long int kmer2_incidence = smooth_kmer_counts(results, file_number, compared_kmer, current_kmer_length, current_gap_position, current_gap_length, background_subtraction);
    long int background;

    if (background_subtraction) {
        if (scale_factor >= 1.0) {
            background = smooth_kmer_counts(results, 0, compared_kmer, current_kmer_length, current_gap_position, current_gap_length, background_subtraction);
            kmer2_incidence -= background * (long int)scale_factor;
        } else {
            background = smooth_kmer_counts(results, 0, compared_kmer, current_kmer_length, current_gap_position, current_gap_length, background_subtraction);
            kmer2_incidence = kmer2_incidence * (long int)(1.0/scale_factor) - background;
        }
    }
    
    // Regular comparison for same length kmers
    if (current_kmer_length == original_kmer_length) {
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

short int Localmax(long int *****results, short int file_number, short int current_kmer_length,
                  short int shortest_kmer, short int longest_kmer_counted, short int current_gap_position,
                  short int current_gap_length, long int current_kmer, double kmer_length_difference_cutoff,
                  short int background_subtraction, double scale_factor)
{
    short int original_kmer_length = current_kmer_length;
    short int count_also_spaced_kmers = 1;
    short int too_long_kmer = longest_kmer_counted + 1;

    // Calculate initial kmer1_incidence
    long int kmer1_incidence;
    if (background_subtraction) {
        if (scale_factor >= 1.0) {
            long int background = smooth_kmer_counts(results, 0, current_kmer, current_kmer_length, current_gap_position, current_gap_length, 1);
            kmer1_incidence = smooth_kmer_counts(results, file_number, current_kmer, current_kmer_length, current_gap_position, current_gap_length, background_subtraction)
                             - background * (long int)scale_factor;
        } else {
            long int background = smooth_kmer_counts(results, 0, current_kmer, current_kmer_length, current_gap_position, current_gap_length, 1);
            kmer1_incidence = smooth_kmer_counts(results, file_number, current_kmer, current_kmer_length, current_gap_position, current_gap_length, background_subtraction)
                             * (long int)(1.0/scale_factor) - background;
        }
    } else {
        kmer1_incidence = smooth_kmer_counts(results, file_number, current_kmer, current_kmer_length, current_gap_position, current_gap_length, background_subtraction);
    }

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
    
    /* Substitution; HAMMING OF 1 */
    for(position=0; position < current_kmer_length; position++, lowbit <<= 2, highbit <<= 2) {
        compared_kmer = lowbit ^ current_kmer;
        if(check_kmer(results, file_number, current_kmer_length, current_gap_position, current_gap_length, compared_kmer, kmer1_incidence, background_subtraction, scale_factor, kmer_length_difference_cutoff, original_kmer_length)) return 0;
        
        compared_kmer = highbit ^ current_kmer;
        if(check_kmer(results, file_number, current_kmer_length, current_gap_position, current_gap_length, compared_kmer, kmer1_incidence, background_subtraction, scale_factor, kmer_length_difference_cutoff, original_kmer_length)) return 0;
        
        compared_kmer = lowbit ^ compared_kmer;
        if(check_kmer(results, file_number, current_kmer_length, current_gap_position, current_gap_length, compared_kmer, kmer1_incidence, background_subtraction, scale_factor, kmer_length_difference_cutoff, original_kmer_length)) return 0;
    }
    
    /* Shift */
    shift = (current_gap_length != 0);
    current_gap_position += shift;
    
    if (current_gap_length == 0 || (count_also_spaced_kmers == 2 - SHIFTED_GAP_POSITIONS && current_gap_position < current_kmer_length)
        || (current_kmer_length % 2 == 1 && count_also_spaced_kmers == 1 && current_gap_position == current_kmer_length / 2 + 1)) {
        
        compared_kmer = (current_kmer >> 2) & lowmask_ULL[current_kmer_length-1];
        if(check_kmer(results, file_number, current_kmer_length, current_gap_position, current_gap_length, compared_kmer, kmer1_incidence, background_subtraction, scale_factor, kmer_length_difference_cutoff, original_kmer_length)) return 0;
        
        lowbit >>= 2; highbit >>= 2;
        compared_kmer = lowbit ^ compared_kmer;
        if(check_kmer(results, file_number, current_kmer_length, current_gap_position, current_gap_length, compared_kmer, kmer1_incidence, background_subtraction, scale_factor, kmer_length_difference_cutoff, original_kmer_length)) return 0;
        
        compared_kmer = highbit ^ compared_kmer;
        if(check_kmer(results, file_number, current_kmer_length, current_gap_position, current_gap_length, compared_kmer, kmer1_incidence, background_subtraction, scale_factor, kmer_length_difference_cutoff, original_kmer_length)) return 0;
        
        compared_kmer = lowbit ^ compared_kmer;
        if(check_kmer(results, file_number, current_kmer_length, current_gap_position, current_gap_length, compared_kmer, kmer1_incidence, background_subtraction, scale_factor, kmer_length_difference_cutoff, original_kmer_length)) return 0;
    }
    
    current_gap_position -= shift;
    current_gap_position -= shift;
    if (current_gap_length == 0 || (count_also_spaced_kmers == 2 - SHIFTED_GAP_POSITIONS && current_gap_position > 0)
        || (current_kmer_length % 2 == 1 && count_also_spaced_kmers == 1 && current_gap_position == current_kmer_length / 2)) {
        
        compared_kmer = (current_kmer << 2) & lowmask_ULL[current_kmer_length-1];
        if(check_kmer(results, file_number, current_kmer_length, current_gap_position, current_gap_length, compared_kmer, kmer1_incidence, background_subtraction, scale_factor, kmer_length_difference_cutoff, original_kmer_length)) return 0;
        
        lowbit = 1; highbit = 2;
        compared_kmer = lowbit ^ compared_kmer;
        if(check_kmer(results, file_number, current_kmer_length, current_gap_position, current_gap_length, compared_kmer, kmer1_incidence, background_subtraction, scale_factor, kmer_length_difference_cutoff, original_kmer_length)) return 0;
        
        compared_kmer = highbit ^ compared_kmer;
        if(check_kmer(results, file_number, current_kmer_length, current_gap_position, current_gap_length, compared_kmer, kmer1_incidence, background_subtraction, scale_factor, kmer_length_difference_cutoff, original_kmer_length)) return 0;
        
        compared_kmer = lowbit ^ compared_kmer;
        if(check_kmer(results, file_number, current_kmer_length, current_gap_position, current_gap_length, compared_kmer, kmer1_incidence, background_subtraction, scale_factor, kmer_length_difference_cutoff, original_kmer_length)) return 0;
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
                if(check_kmer(results, file_number, current_kmer_length, current_gap_position, current_gap_length, compared_kmer, kmer1_incidence, background_subtraction, scale_factor, kmer_length_difference_cutoff, original_kmer_length)) return 0;
                
                compared_kmer = lowbit ^ compared_kmer;
                if(check_kmer(results, file_number, current_kmer_length, current_gap_position, current_gap_length, compared_kmer, kmer1_incidence, background_subtraction, scale_factor, kmer_length_difference_cutoff, original_kmer_length)) return 0;
                
                compared_kmer = highbit ^ compared_kmer;
                if(check_kmer(results, file_number, current_kmer_length, current_gap_position, current_gap_length, compared_kmer, kmer1_incidence, background_subtraction, scale_factor, kmer_length_difference_cutoff, original_kmer_length)) return 0;
                
                compared_kmer = lowbit ^ compared_kmer;
                if(check_kmer(results, file_number, current_kmer_length, current_gap_position, current_gap_length, compared_kmer, kmer1_incidence, background_subtraction, scale_factor, kmer_length_difference_cutoff, original_kmer_length)) return 0;
                
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
                if(check_kmer(results, file_number, current_kmer_length, current_gap_position, current_gap_length, compared_kmer, kmer1_incidence, background_subtraction, scale_factor, kmer_length_difference_cutoff, original_kmer_length)) return 0;
                
                compared_kmer = lowbit ^ compared_kmer;
                if(check_kmer(results, file_number, current_kmer_length, current_gap_position, current_gap_length, compared_kmer, kmer1_incidence, background_subtraction, scale_factor, kmer_length_difference_cutoff, original_kmer_length)) return 0;
                
                compared_kmer = highbit ^ compared_kmer;
                if(check_kmer(results, file_number, current_kmer_length, current_gap_position, current_gap_length, compared_kmer, kmer1_incidence, background_subtraction, scale_factor, kmer_length_difference_cutoff, original_kmer_length)) return 0;
                
                compared_kmer = lowbit ^ compared_kmer;
                if(check_kmer(results, file_number, current_kmer_length, current_gap_position, current_gap_length, compared_kmer, kmer1_incidence, background_subtraction, scale_factor, kmer_length_difference_cutoff, original_kmer_length)) return 0;
                
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
                            
                            if(check_kmer(results, file_number, current_kmer_length, current_gap_position, current_gap_length, compared_kmer, kmer1_incidence, background_subtraction, scale_factor, kmer_length_difference_cutoff, original_kmer_length)) return 0;
                            compared_kmer = lowbit ^ compared_kmer;
                            if(check_kmer(results, file_number, current_kmer_length, current_gap_position, current_gap_length, compared_kmer, kmer1_incidence, background_subtraction, scale_factor, kmer_length_difference_cutoff, original_kmer_length)) return 0;
                            compared_kmer = highbit ^ compared_kmer;
                            if(check_kmer(results, file_number, current_kmer_length, current_gap_position, current_gap_length, compared_kmer, kmer1_incidence, background_subtraction, scale_factor, kmer_length_difference_cutoff, original_kmer_length)) return 0;
                            compared_kmer = lowbit ^ compared_kmer;
                            if(check_kmer(results, file_number, current_kmer_length, current_gap_position, current_gap_length, compared_kmer, kmer1_incidence, background_subtraction, scale_factor, kmer_length_difference_cutoff, original_kmer_length)) return 0;
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
                            if(check_kmer(results, file_number, current_kmer_length, current_gap_position, current_gap_length, compared_kmer, kmer1_incidence, background_subtraction, scale_factor, kmer_length_difference_cutoff, original_kmer_length)) return 0;
                            compared_kmer = lowbit ^ compared_kmer;
                            if(check_kmer(results, file_number, current_kmer_length, current_gap_position, current_gap_length, compared_kmer, kmer1_incidence, background_subtraction, scale_factor, kmer_length_difference_cutoff, original_kmer_length)) return 0;
                            compared_kmer = highbit ^ compared_kmer;
                            if(check_kmer(results, file_number, current_kmer_length, current_gap_position, current_gap_length, compared_kmer, kmer1_incidence, background_subtraction, scale_factor, kmer_length_difference_cutoff, original_kmer_length)) return 0;
                            compared_kmer = lowbit ^ compared_kmer;
                            if(check_kmer(results, file_number, current_kmer_length, current_gap_position, current_gap_length, compared_kmer, kmer1_incidence, background_subtraction, scale_factor, kmer_length_difference_cutoff, original_kmer_length)) return 0;
                            
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
                            if(check_kmer(results, file_number, current_kmer_length, current_gap_position, current_gap_length, compared_kmer, kmer1_incidence, background_subtraction, scale_factor, kmer_length_difference_cutoff, original_kmer_length)) return 0;
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
                        if(check_kmer(results, file_number, current_kmer_length, current_gap_position, current_gap_length, compared_kmer, kmer1_incidence, background_subtraction, scale_factor, kmer_length_difference_cutoff, original_kmer_length)) return 0;
                    }
                }
                
                /* Right part */
                if (current_gap_position != 1 && right == 1) {
                    if (current_gap_position > 0) current_gap_position--;
                    compared_kmer = (current_kmer & lowmask_ULL[current_kmer_length-1]);
                    if(check_kmer(results, file_number, current_kmer_length, current_gap_position, current_gap_length, compared_kmer, kmer1_incidence, background_subtraction, scale_factor, kmer_length_difference_cutoff, original_kmer_length)) return 0;
                }
                
                current_gap_position = true_gap_position;
                
                /* Shorter with longer gap */
                if(count_also_spaced_kmers != 0 && true_gap_position != 0 && current_gap_length < MAX_GAP_LENGTH) {
                    current_gap_length++;
                    if (current_gap_position < current_kmer_length && left == 1) {
                        compared_kmer = ((current_kmer >> 2) & highmask_ULL[current_kmer_length - current_gap_position]) | ((current_kmer) & (lowmask_ULL[current_kmer_length-current_gap_position-1]));
                        if(check_kmer(results, file_number, current_kmer_length, current_gap_position, current_gap_length, compared_kmer, kmer1_incidence, background_subtraction, scale_factor, kmer_length_difference_cutoff, original_kmer_length)) return 0;
                    }
                    
                    if(current_gap_position > 1 && right == 1) {
                        current_gap_position--;
                        compared_kmer = ((current_kmer >> 2) & highmask_ULL[current_kmer_length - current_gap_position]) | ((current_kmer) & (lowmask_ULL[current_kmer_length-current_gap_position-1]));
                        if(check_kmer(results, file_number, current_kmer_length, current_gap_position, current_gap_length, compared_kmer, kmer1_incidence, background_subtraction, scale_factor, kmer_length_difference_cutoff, original_kmer_length)) return 0;
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
                    if(check_kmer(results, file_number, current_kmer_length, current_gap_position, current_gap_length, compared_kmer, kmer1_incidence, background_subtraction, scale_factor, kmer_length_difference_cutoff, original_kmer_length)) return 0;
                }
                else if(current_kmer_length-current_gap_position == 0) {
                    current_gap_length = 0;
                    current_gap_position = 0;
                    compared_kmer = (current_kmer >> 2);
                    if(check_kmer(results, file_number, current_kmer_length, current_gap_position, current_gap_length, compared_kmer, kmer1_incidence, background_subtraction, scale_factor, kmer_length_difference_cutoff, original_kmer_length)) return 0;
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
                        
                        if(check_kmer(results, file_number, current_kmer_length, current_gap_position, current_gap_length, compared_kmer, kmer1_incidence, background_subtraction, scale_factor, kmer_length_difference_cutoff, original_kmer_length)) return 0;
                        compared_kmer = lowbit ^ compared_kmer;
                        if(check_kmer(results, file_number, current_kmer_length, current_gap_position, current_gap_length, compared_kmer, kmer1_incidence, background_subtraction, scale_factor, kmer_length_difference_cutoff, original_kmer_length)) return 0;
                        compared_kmer = highbit ^ compared_kmer;
                        if(check_kmer(results, file_number, current_kmer_length, current_gap_position, current_gap_length, compared_kmer, kmer1_incidence, background_subtraction, scale_factor, kmer_length_difference_cutoff, original_kmer_length)) return 0;
                        compared_kmer = lowbit ^ compared_kmer;
                        if(check_kmer(results, file_number, current_kmer_length, current_gap_position, current_gap_length, compared_kmer, kmer1_incidence, background_subtraction, scale_factor, kmer_length_difference_cutoff, original_kmer_length)) return 0;
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
                        if(count_also_spaced_kmers != 0 && true_gap_length >= 1) {
                            current_gap_length--;
                            lowbit = 1; highbit = 2;
                            lowbit <<= ((current_kmer_length - current_gap_position-1)*2); highbit <<= ((current_kmer_length - current_gap_position-1)*2);
                            compared_kmer = ((current_kmer << 2) & highmask_ULL[current_kmer_length - current_gap_position]) | ((current_kmer) & (lowmask_ULL[current_kmer_length-current_gap_position-1]));
                            if (current_gap_length == 0) current_gap_position = 0;
                            
                            for(first_half = 0; first_half < 2; first_half++) {
                                if (current_gap_position < current_kmer_length && (current_gap_position == 0 || Centergap (count_also_spaced_kmers, current_kmer_length, current_gap_position))) {
                                    if(check_kmer(results, file_number, current_kmer_length, current_gap_position, current_gap_length, compared_kmer, kmer1_incidence, background_subtraction, scale_factor, kmer_length_difference_cutoff, original_kmer_length)) return 0;
                                    compared_kmer = lowbit ^ compared_kmer;
                                    if(check_kmer(results, file_number, current_kmer_length, current_gap_position, current_gap_length, compared_kmer, kmer1_incidence, background_subtraction, scale_factor, kmer_length_difference_cutoff, original_kmer_length)) return 0;
                                    compared_kmer = highbit ^ compared_kmer;
                                    if(check_kmer(results, file_number, current_kmer_length, current_gap_position, current_gap_length, compared_kmer, kmer1_incidence, background_subtraction, scale_factor, kmer_length_difference_cutoff, original_kmer_length)) return 0;
                                    compared_kmer = lowbit ^ compared_kmer;
                                    if(check_kmer(results, file_number, current_kmer_length, current_gap_position, current_gap_length, compared_kmer, kmer1_incidence, background_subtraction, scale_factor, kmer_length_difference_cutoff, original_kmer_length)) return 0;
                                }
                                
                                current_gap_position++;
                                if (current_gap_position > current_kmer_length || current_gap_position == 1) break;
                                compared_kmer = ((current_kmer << 2) & highmask_ULL[current_kmer_length - current_gap_position]) | ((current_kmer) & (lowmask_ULL[current_kmer_length-current_gap_position-1]));
                            }
                            current_gap_length++;
                        }
                    }
                    current_gap_position = true_gap_position;
                    return 1;
                }
        
/* SUBROUTINE THAT DETERMINES IF A GAPPED KMER IS A LOCAL MAXIMUM WITHIN HUDDINGE DISTANCE OF 1 */
/* SEE NITTA ET AL. eLIFE 2015 Methods and Supplementary Figure 1 for algorithm description     */
/* https://doi.org/10.7554/eLife.04837.004                                                      */
/* Bases encoded as bits, A = 00, C = 01, G = 10, T = 11                                        */

short int oldLocalmax(long int *****results, short int file_number, short int current_kmer_length, short int shortest_kmer, short int longest_kmer_counted, short int current_gap_position, short int current_gap_length, long int current_kmer, double kmer_length_difference_cutoff, short int noeffect, double another_noeffect)
{

    short int count_also_spaced_kmers = 1;
    short int too_long_kmer = longest_kmer_counted + 1;

    long int kmer1_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][current_kmer];
    long int kmer2_incidence;
    long int compared_kmer = current_kmer;
    signed short int position;
    short int counter;
    short int first_half;
    long int lowbit = 1;  // "01" used to toggle low bit of base encoding (between A-C or G-T), shifted to correct base position in loops
    long int highbit = 2; // "10" used to toggle high bit of base encoding (between A-G or C-T), shifted to correct base position in loops
    short int shift;
    short int true_gap_position = current_gap_position;
    short int true_gap_length = current_gap_length;
    short int start;
    short int end;
    short int left = 0;
    short int right = 1;
    signed short int position2;
    
    /* Substitution; HAMMING OF 1, RETURNS 0 IF ANY KMER WITHIN HAMMING OF 1 HAS HIGHER COUNT */
    for(position=0; position < current_kmer_length; position++, lowbit <<= 2, highbit <<= 2)
    {
        compared_kmer = lowbit ^ current_kmer;  // SHIFTS between A-C or G-T of original kmer at position
        kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
        if(kmer2_incidence > kmer1_incidence) return (0);  // IF SUBSTITUTED KMER HAS HIGHER COUNTS, RETURNS 0 TO INDICATE THAT THE QUERY KMER IS NOT LOCAL MAXIMUM
        compared_kmer = highbit ^ current_kmer; // SHIFTS between A-G or C-T of original kmer at position
        kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
        if(kmer2_incidence > kmer1_incidence) return (0);
        compared_kmer = lowbit ^ compared_kmer; // SHIFTS between A-C or G-T of previously shifted kmer to cover last of possible 3 substitutions (total effect is A-T, C-G)
        kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
        if(kmer2_incidence > kmer1_incidence) return (0);
    }
    
    /* Shift; FULL SHIFT BY ONE */
    shift = (current_gap_length != 0);
    current_gap_position += shift;
    
    /* ONLY LOOK AT FULL SHIFT FOR UNGAPPED KMERS, OR FOR GAPPED KMERS IF GAP CAN SHIFT ALSO (ALL GAP POSITIONS HAVE BEEN COUNTED) */
    // original :: if (current_gap_length == 0 || (count_also_spaced_kmers == 2 && current_gap_position < current_kmer_length) || (current_kmer_length % 2 == 1 && count_also_spaced_kmers == 1 && current_gap_position == current_kmer_length / 2 + 1))
    if (current_gap_length == 0 || (count_also_spaced_kmers == 2 - SHIFTED_GAP_POSITIONS && current_gap_position < current_kmer_length) || (current_kmer_length % 2 == 1 && count_also_spaced_kmers == 1 && current_gap_position == current_kmer_length / 2 + 1)) // looks at
    {
    /* COMPARED FULL SHIFT RIGHT (same shift in gap if any), RETURNS 0 IF ANY KMER WITHIN HAMMING OF 1 HAS HIGHER COUNT  */
    compared_kmer = (current_kmer >> 2) & lowmask_ULL[current_kmer_length-1];
    kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
    if(kmer2_incidence > kmer1_incidence) return (0);
    lowbit >>= 2; highbit >>= 2;
    compared_kmer = lowbit ^ compared_kmer;
    kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
    if(kmer2_incidence > kmer1_incidence) return (0);
    compared_kmer = highbit ^ compared_kmer;
    kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
    if(kmer2_incidence > kmer1_incidence) return (0);
    compared_kmer = lowbit ^ compared_kmer;
    kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
    if(kmer2_incidence > kmer1_incidence) return (0);
    }
    
    current_gap_position -= shift;
    current_gap_position -= shift;
    // if (current_gap_length == 0 || (count_also_spaced_kmers == 2 && current_gap_position > 0) || (current_kmer_length % 2 == 1 && count_also_spaced_kmers == 1 && current_gap_position == current_kmer_length / 2))
    if (current_gap_length == 0 || (count_also_spaced_kmers == 2 - SHIFTED_GAP_POSITIONS && current_gap_position > 0) || (current_kmer_length % 2 == 1 && count_also_spaced_kmers == 1 && current_gap_position == current_kmer_length / 2))
    {
    /* COMPARED FULL SHIFT LEFT (same shift in gap if any) */
    compared_kmer = (current_kmer << 2) & lowmask_ULL[current_kmer_length-1];
    kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
    if(kmer2_incidence > kmer1_incidence) return (0);
    lowbit = 1; highbit = 2;
    compared_kmer = lowbit ^ compared_kmer;
    kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
    if(kmer2_incidence > kmer1_incidence) return (0);
    compared_kmer = highbit ^ compared_kmer;
    kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
    if(kmer2_incidence > kmer1_incidence) return (0);
    compared_kmer = lowbit ^ compared_kmer;
    kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
    if(kmer2_incidence > kmer1_incidence) return (0);
    }
    
    current_gap_position += shift;
    lowbit = 1; highbit = 2;
    if (count_also_spaced_kmers != 0)
    {
    if(current_gap_position == 0) current_gap_position = current_kmer_length / 2;
    
    /* Longer Gap; COMPARE TO KMER WITH LONGER GAP */
    current_gap_length++;
    if (current_gap_length < Nlength - current_kmer_length && current_gap_length < MAX_GAP_LENGTH && true_gap_length != 0)
    {
    compared_kmer = (current_kmer & highmask_ULL[current_kmer_length-current_gap_position]) | ((current_kmer << 2) & (lowmask_ULL[current_kmer_length-current_gap_position-1]));
    
    /* LOOP TO ANALYZE SHIFT OF EITHER HALF, STARTS WITH SECOND HALF */
    for(first_half = 0; first_half < 2; first_half++)
    {
    kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
    if(kmer2_incidence > kmer1_incidence) return (0);
    compared_kmer = lowbit ^ compared_kmer;
    kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
    if(kmer2_incidence > kmer1_incidence) return (0);
    compared_kmer = highbit ^ compared_kmer;
    kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
    if(kmer2_incidence > kmer1_incidence) return (0);
    compared_kmer = lowbit ^ compared_kmer;
    kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
    if(kmer2_incidence > kmer1_incidence) return (0);
    
    /* SWITCHES TO FIRST HALF */
    compared_kmer = ((current_kmer >> 2) & highmask_ULL[current_kmer_length-current_gap_position]) | ((current_kmer) & (lowmask_ULL[current_kmer_length-current_gap_position-1]));
    lowbit <<= ((current_kmer_length - 1)*2); highbit <<= ((current_kmer_length - 1)*2);
    }
    }
        
    /* Shorter Gap; COMPARE TO KMER WITH SHORTER GAP */
    current_gap_length--;
    current_gap_length--;

    lowbit = 1; highbit = 2;
    lowbit <<= ((current_kmer_length - true_gap_position - 1)*2); highbit <<= ((current_kmer_length - true_gap_position - 1)*2);
    if (current_gap_length == 0) current_gap_position = 0;
    if (current_gap_length >= 0)
    {
        compared_kmer = (current_kmer & highmask_ULL[current_kmer_length - true_gap_position]) | ((current_kmer >> 2) & (lowmask_ULL[current_kmer_length-true_gap_position-1]));
        /* LOOP TO ANALYZE SHIFT OF EITHER HALF, STARTS WITH SECOND HALF */
        for(first_half = 0; first_half < 2; first_half++)
        {
            kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
            if(kmer2_incidence > kmer1_incidence) return (0);
            compared_kmer = lowbit ^ compared_kmer;
            kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
            if(kmer2_incidence > kmer1_incidence) return (0);
            compared_kmer = highbit ^ compared_kmer;
            kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
            if(kmer2_incidence > kmer1_incidence) return (0);
            compared_kmer = lowbit ^ compared_kmer;
            kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
            if(kmer2_incidence > kmer1_incidence) return (0);
            
            /* SWITCHES TO FIRST HALF */
            compared_kmer = lowmask_ULL[current_kmer_length-1] & ((current_kmer << 2) & highmask_ULL[current_kmer_length - true_gap_position]) | ((current_kmer) & (lowmask_ULL[current_kmer_length-true_gap_position-1]));
            lowbit <<= 2; highbit <<= 2;
        }
    }
        
        /* COMPARES DIFFERENT GAP POSITIONS (SHIFTED BY ONE) */
        current_gap_length = true_gap_length;
        if ((count_also_spaced_kmers == 2 || (count_also_spaced_kmers == 1 && current_kmer_length % 2 == 1)) && true_gap_length > 0)
        {
        current_gap_position = true_gap_position + 1;
        lowbit = 1; highbit = 2;
        lowbit <<= ((current_kmer_length - current_gap_position)*2); highbit <<= ((current_kmer_length - current_gap_position)*2);
        compared_kmer = (current_kmer & highmask_ULL[current_kmer_length - current_gap_position]) | ((current_kmer) & (lowmask_ULL[current_kmer_length-current_gap_position-1]));

        /* LOOP TO ANALYZE SHIFT TO EITHER DIRECTION, STARTS WITH RIGHT BY ONE */
        for(first_half = 0; first_half < 2; first_half++)
        {
        if (current_gap_position < current_kmer_length && current_gap_position > 0 && (count_also_spaced_kmers == 2 ||
        (count_also_spaced_kmers == 1 && ((current_gap_position == current_kmer_length / 2) || current_gap_position == current_kmer_length / 2 + 1))))
        {
        kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
        if(kmer2_incidence > kmer1_incidence) return (0);
        compared_kmer = lowbit ^ compared_kmer;
        kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
        if(kmer2_incidence > kmer1_incidence) return (0);
        compared_kmer = highbit ^ compared_kmer;
        kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
        if(kmer2_incidence > kmer1_incidence) return (0);
        compared_kmer = lowbit ^ compared_kmer;
        kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
        if(kmer2_incidence > kmer1_incidence) return (0);
        }
            
        /* SWITCHES LEFT BY ONE */
        current_gap_position--;
        current_gap_position--;
        compared_kmer = ((current_kmer) & highmask_ULL[current_kmer_length - current_gap_position]) | ((current_kmer) & (lowmask_ULL[current_kmer_length-current_gap_position-1]));
        lowbit <<= 2; highbit <<= 2;
        }
        }
        
        start = 1;
        end = current_kmer_length;
        
        /* COMPARES UNGAPPED KMER TO ALL SINGLE GAPS IN ALL POSITIONS */
        if (count_also_spaced_kmers != 0 && true_gap_length == 0)
        {
            if(count_also_spaced_kmers == 1)
            {
                start = current_kmer_length / 2;
                end = start + 1 + (current_kmer_length % 2);
            }
        for(current_gap_position = start, current_gap_length = 1; current_gap_position < end; current_gap_position++)
        {
        compared_kmer = (current_kmer & highmask_ULL[current_kmer_length - current_gap_position]) | ((current_kmer << 2) & (lowmask_ULL[current_kmer_length-current_gap_position-1]));
            /* LOOP TO ANALYZE SHIFT TO EITHER DIRECTION, STARTS WITH BEGINNING */
            for(lowbit = 1, highbit = 2, first_half = 0; first_half < 2; first_half++)
            {
            kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
            if(kmer2_incidence > kmer1_incidence) return (0);
            compared_kmer = lowbit ^ compared_kmer;
            kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
            if(kmer2_incidence > kmer1_incidence) return (0);
            compared_kmer = highbit ^ compared_kmer;
            kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
            if(kmer2_incidence > kmer1_incidence) return (0);
            compared_kmer = lowbit ^ compared_kmer;
            kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
            if(kmer2_incidence > kmer1_incidence) return (0);
           
        /* END */
        compared_kmer = ((current_kmer >> 2) & highmask_ULL[current_kmer_length - current_gap_position]) | ((current_kmer) & (lowmask_ULL[current_kmer_length-current_gap_position-1]));
            lowbit <<= ((current_kmer_length - 1)*2); highbit <<= ((current_kmer_length - 1)*2);
        }
        }
        }
        
    }
/* Shorter; COMPARES KMER WITH ONE SHORTER */
current_gap_position = true_gap_position;
current_gap_length = true_gap_length;
current_kmer_length--;
end = current_kmer_length;
if (current_kmer_length >= shortest_kmer)
{
/* IF NO GAP, INSERTS GAP AT ALL ALLOWED POSITIONS */
if(current_gap_length == 0)
{
    if (count_also_spaced_kmers != 0)
    {
        if(count_also_spaced_kmers == 1)
        {
            start = current_kmer_length / 2;
            end = start + 1 + (current_kmer_length % 2);
        }
        for(current_gap_position = start, current_gap_length = 1; current_gap_position < end; current_gap_position++)
        {
            compared_kmer = ((current_kmer >> 2) & highmask_ULL[current_kmer_length - current_gap_position]) | ((current_kmer) & (lowmask_ULL[current_kmer_length-current_gap_position-1]));
                kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
                if(kmer2_incidence  * kmer_length_difference_cutoff > kmer1_incidence) return (0);
        }
    }
}

current_gap_position = true_gap_position;
current_gap_length = true_gap_length;
    
if (count_also_spaced_kmers != 1) {left = 1; right = 1;}
if (current_gap_position == current_kmer_length / 2 && current_kmer_length % 2 == 0 && count_also_spaced_kmers == 1) {left = 1; right = 0;}
if (current_kmer_length % 2 == 1 && count_also_spaced_kmers == 1) left = 1;
    
/* LEFT PART */
if (current_gap_position < current_kmer_length)
{
    if (left == 1 || true_gap_length == 0)
    {
compared_kmer = (current_kmer >> 2);
kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
if(kmer2_incidence * kmer_length_difference_cutoff > kmer1_incidence) return (0);
    }
    }
/* RIGHT PART */
if (current_gap_position != 1 && right == 1)
{
if (current_gap_position > 0) current_gap_position--;
    compared_kmer = (current_kmer & lowmask_ULL[current_kmer_length-1]);
kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
if(kmer2_incidence * kmer_length_difference_cutoff > kmer1_incidence) return (0);
}
current_gap_position = true_gap_position;
/* Shorter with Longer Gap; LONGER GAP */
if(count_also_spaced_kmers != 0 && true_gap_position != 0 && current_gap_length < MAX_GAP_LENGTH)
{
current_gap_length++;
if (current_gap_position < current_kmer_length && left == 1)
{
/* RIGHT BASE GAPPED */
compared_kmer = ((current_kmer >> 2) & highmask_ULL[current_kmer_length - current_gap_position]) | ((current_kmer) & (lowmask_ULL[current_kmer_length-current_gap_position-1]));
    /* printf("DEBUG Localmax: k=%d, gap_pos=%d, gap_len=%d, indices: high=%d, low=%d\n",
           current_kmer_length, current_gap_position, current_gap_length,
           current_kmer_length - current_gap_position,
           current_kmer_length-current_gap_position-1); fflush(stdout); */
kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
if(kmer2_incidence * kmer_length_difference_cutoff > kmer1_incidence) return (0);
}
/* LEFT BASE GAPPED */
if(current_gap_position > 1 && right == 1)
{
current_gap_position--;
compared_kmer = ((current_kmer >> 2) & highmask_ULL[current_kmer_length - current_gap_position]) | ((current_kmer) & (lowmask_ULL[current_kmer_length-current_gap_position-1]));
kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
if(kmer2_incidence * kmer_length_difference_cutoff > kmer1_incidence) return (0);
}
current_gap_length--;
current_gap_position++;
}
    /* COMPARES HANGING SINGLE BASE TO UNGAPPED KMER */
    current_gap_position = true_gap_position;
    current_gap_length = true_gap_length;
    if(current_gap_position == 1)
    {
    current_gap_length = 0;
    current_gap_position = 0;
    compared_kmer = (current_kmer & lowmask_ULL[current_kmer_length-1]);
    kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
    if(kmer2_incidence * kmer_length_difference_cutoff > kmer1_incidence) return (0);
    }
    else if(current_kmer_length-current_gap_position == 0)
    {
    current_gap_length = 0;
    current_gap_position = 0;
    compared_kmer = (current_kmer >> 2);
    kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
    if(kmer2_incidence * kmer_length_difference_cutoff > kmer1_incidence) return (0);
    }
}
current_gap_position = true_gap_position;
current_gap_length = true_gap_length;
    
    /* Longer; COMPARES KMER WITH ONE LONGER */
    current_kmer_length++;
    current_kmer_length++;
    if (current_kmer_length < too_long_kmer)
    {
        //if (current_kmer == 0xE34) {  printf("\nDEBUG: TGAnTCA (gap=%d) comparing with longer\n", current_gap_position);}
        compared_kmer = (current_kmer << 2);
        /* LOOP TO ANALYZE SHIFT TO EITHER DIRECTION, STARTS WITH BEGINNING */
        for(lowbit = 1, highbit = 2, first_half = 0; first_half < 2; first_half++)
        {
            if(count_also_spaced_kmers != 1 || true_gap_length == 0 ||
               (first_half == 1 || (count_also_spaced_kmers == 1 && current_gap_position == current_kmer_length / 2)))
            {
                kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
                    // if(current_kmer == 0xE34) {
                    // printf("\nDEBUG compare: kmer1 ");
                    // Kmerprint(compared_kmer, current_kmer_length, current_gap_position, current_gap_length);
                    // printf("\tcount=%ld kmer2_count=%ld threshold=%f length=%d gap_pos=%d\n", kmer1_incidence, kmer2_incidence, kmer1_incidence * kmer_length_difference_cutoff, current_kmer_length, current_gap_position);}
                if(kmer2_incidence > kmer1_incidence * kmer_length_difference_cutoff) return (0);
                compared_kmer = lowbit ^ compared_kmer;
                kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
                
                if(kmer2_incidence > kmer1_incidence * kmer_length_difference_cutoff) return (0);
                compared_kmer = highbit ^ compared_kmer;
                kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
                
                if(kmer2_incidence > kmer1_incidence * kmer_length_difference_cutoff) return (0);
                compared_kmer = lowbit ^ compared_kmer;
                kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
                
                if(kmer2_incidence > kmer1_incidence * kmer_length_difference_cutoff) return (0);
            }
            /* END */
            compared_kmer = current_kmer;
            lowbit <<= ((current_kmer_length - 1)*2); highbit <<= ((current_kmer_length - 1)*2);
            if (true_gap_length == 0) continue;
            current_gap_position++;
            if (current_gap_position > current_kmer_length || (count_also_spaced_kmers == 1 && current_gap_position != (current_kmer_length / 2 + (current_kmer_length % 2)))) break;
        }
        current_gap_position = true_gap_position;
        current_gap_length = true_gap_length;
    

    
/* Longer with Shorter Gap; COMPARES TO ONE SHORTER GAP LENGTH */
if(count_also_spaced_kmers != 0 && true_gap_length >= 1)
{
current_gap_length--;
lowbit = 1; highbit = 2;
lowbit <<= ((current_kmer_length - current_gap_position-1)*2); highbit <<= ((current_kmer_length - current_gap_position-1)*2);
compared_kmer = ((current_kmer << 2) & highmask_ULL[current_kmer_length - current_gap_position]) | ((current_kmer) & (lowmask_ULL[current_kmer_length-current_gap_position-1]));
/* NO GAP LEFT */
if (current_gap_length == 0) current_gap_position = 0;
    
    /* LOOP TO ANALYZE ADDED BASE ON EITHER SIDE OF GAP, STARTS WITH RIGHT */
    for(first_half = 0; first_half < 2; first_half++)
    {
        if (current_gap_position < current_kmer_length && (current_gap_position == 0 || Centergap (count_also_spaced_kmers, current_kmer_length, current_gap_position)))
        {
            kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
            if(kmer2_incidence > kmer1_incidence * kmer_length_difference_cutoff) return (0);
            compared_kmer = lowbit ^ compared_kmer;
            kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
            if(kmer2_incidence > kmer1_incidence * kmer_length_difference_cutoff) return (0);
            compared_kmer = highbit ^ compared_kmer;
            kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
            if(kmer2_incidence > kmer1_incidence * kmer_length_difference_cutoff) return (0);
            compared_kmer = lowbit ^ compared_kmer;
            kmer2_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][compared_kmer];
            if(kmer2_incidence > kmer1_incidence * kmer_length_difference_cutoff) return (0);
        }
        
        /* SWITCHES ADDED BASE TO LEFT */
        current_gap_position++;
        if (current_gap_position > current_kmer_length || current_gap_position == 1) break;
        compared_kmer = ((current_kmer << 2) & highmask_ULL[current_kmer_length - current_gap_position]) | ((current_kmer) & (lowmask_ULL[current_kmer_length-current_gap_position-1]));
    }
current_gap_length++;
}
}
current_gap_position = true_gap_position;
return (1);
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

// Function declarations for externally defined functions
// Function to encode a DNA character to 2-bit representation.
uint64_t encode_base(char base) {
    switch (base) {
        case 'A': return 0x0ULL; // 00
        case 'C': return 0x1ULL; // 01
        case 'G': return 0x2ULL; // 10
        case 'T': return 0x3ULL; // 11
        default: {/*printf("\nError in DNA sequence");*/return UINT64_MAX;} // Return max if invalid character
    }
}

// Function to convert a kmer index back to a DNA string.
void decode_kmer(uint64_t kmer, size_t kmer_length, size_t gap_position, size_t gap_length, char *out) {
    int Ncount = gap_length;
    int kmer_position = kmer_length;
    int i = kmer_position + gap_length - 1;

    for (; i >= 0; i--) {
        if (kmer_position == gap_position && Ncount > 0) {
            out[i] = 'n';
            Ncount--;
        }
        else {
            switch (kmer & 0x3) {
                case 0x0: out[i] = 'A'; break;
                case 0x1: out[i] = 'C'; break;
                case 0x2: out[i] = 'G'; break;
                case 0x3: out[i] = 'T'; break;
            }
            kmer >>= 2;
            kmer_position--;
        }
    }
    out[kmer_length+gap_length] = '\0'; // Null-terminate the string.
}

void print_bitstream_as_dna(const uint64_t *bitstream, size_t block_index, size_t dna_length) {
    printf("Block %zu: ", block_index);
    for (int i = 31; i >= 0; i--) {
        uint64_t mask = 0x3ULL << (2 * i);
        char base = ((32 - i) < dna_length) ?
                    "ACGT"[(bitstream[block_index] & mask) >> (2 * i)] :
                    'a';
        printf("%c", base);
    }
    printf("\n");
}


short int encode_dna_sequence(const char *dna_sequence, uint64_t *bitstream, size_t *bitstream_length) {
    size_t dna_length = strlen(dna_sequence);
    *bitstream_length = (2 * dna_length + 63) / 64;
    memset(bitstream, 0, (*bitstream_length) * sizeof(uint64_t));
    
    // printf("DEBUG: Encoding sequence of length %zu\n", dna_length);
    
    // Initial encoding
    size_t block;
    for (size_t i = 0; i < dna_length; i++) {
        uint64_t encoded_base = encode_base(dna_sequence[i]);
        if (encoded_base == UINT64_MAX) return (0); // DEFECTIVE SEQUENCE
        block = (*bitstream_length) - 1 - i * 2 / 64;
        size_t offset = 62 - (i * 2 % 64);
        bitstream[block] |= (encoded_base << offset);
    }
    bitstream[(*bitstream_length)] = 0;  // add 0 block after so that shifting last block will not bring in garbage seq

    // printf("\nAfter initial encoding:\n");
    // for (size_t i = 0; i < *bitstream_length; i++) {
        // print_bitstream_as_dna(bitstream, i, dna_length);
    // }

    // Handle short sequences
    if (dna_length <= 32) {
        bitstream[0] = bitstream[0] >> (64 - 2 * dna_length);
        // printf("\nAfter short sequence adjustment:\n");
        // print_bitstream_as_dna(bitstream, 0, dna_length);
    }
    
    // Handle longer sequences
    if (dna_length > 32 && dna_length % 32 != 0) {
        size_t left_shift = (dna_length * 2) % 64;
        size_t right_shift = 64 - left_shift;
        
        /*printf("\nBefore shifting blocks:\n");
        for (size_t i = 0; i < *bitstream_length; i++) {
            print_bitstream_as_dna(bitstream, i, dna_length);
            }*/
        
        uint64_t low_bits;
        //Shift blocks maintaining proper bit transitions
        for (size_t block = 0; block < (*bitstream_length); block++) {
            // Get high bits from next block that need to be brought into this block
            uint64_t high_bits = bitstream[block + 1] << left_shift;
            low_bits = (bitstream[block] >> right_shift);
            // Shift current block right and OR in the high bits from next block
            bitstream[block] = low_bits | high_bits;
        }
        
        /*printf("\nAfter shifting blocks:\n");
        for (size_t i = 0; i < *bitstream_length; i++) {
            print_bitstream_as_dna(bitstream, i, dna_length);
        }*/
        
    }
    return(1);
}

short int reverse_complement_bitstream(uint64_t *bitstream, uint64_t *temp_bitstream,
                                     size_t sequence_length) {
    
    size_t bitstream_length = (2 * sequence_length + 63) / 64;
    
    // Copy input to temp_bitstream first
    for(int i = 0; i < bitstream_length; i++) temp_bitstream[i] = bitstream[i];
 
    // printf("\nInput:\n"); for (size_t i = 0; i < bitstream_length; i++) {print_bitstream_as_dna(bitstream, i, sequence_length);}
    
    // Clear original bitstream for reverse complement
    memset(bitstream, 0, bitstream_length * sizeof(uint64_t));
    
    // Process each base pair (2 bits)
    for (size_t i = 0; i < sequence_length; i++) {
        // Find block and position in input (temp_bitstream)
        size_t in_block = (sequence_length - 1 - i) * 2 / 64;
        size_t in_offset = ((sequence_length - 1 - i) * 2) % 64;
        
        // Get 2 bits from temp_bitstream
        uint64_t base = (temp_bitstream[in_block] >> in_offset) & 0x3;
        
        // Complement the base (00->11, 01->10, 10->01, 11->00)
        base = (~base) & 0x3;
        
        // Calculate output position in original bitstream
        size_t out_block = i * 2 / 64;
        size_t out_offset = (i * 2) % 64;
        
        // Write to original bitstream
        bitstream[out_block] |= (base << out_offset);
    }
    
    //printf("\nOutput:\n"); for (size_t i = 0; i < bitstream_length; i++) {print_bitstream_as_dna(bitstream, i, sequence_length);}
    
    return 1;
}

void count_kmers(uint64_t *results, uint64_t *bitstream, size_t sequence_length, size_t kmer_length, size_t gap_position, size_t gap_length) {
    uint64_t main_reg = 0;
    uint64_t source_reg = 0;
    const int shift_amount = 2;
    size_t total_kmers = sequence_length - kmer_length + 1 - gap_length;
    size_t counted_kmers = 0;
    size_t incoming_bits;
    size_t last_source_block_moved = 0;
    uint64_t kmer;
    size_t current_block = 0;
    // char* kmer_str = (char*)malloc(kmer_length + MAX_GAP_LENGTH + 1);
    // kmer_str[kmer_length] = '\0';
    
    // printf("DEBUG: Starting kmer counting, expecting %zu kmers\n", total_kmers);
    
    size_t bitstream_length = (sequence_length * 2 + 63) / 64;

    uint64_t *current_bitstream = bitstream;

    // Load initial blocks
    if (bitstream_length > 0) {
        main_reg = *current_bitstream++;
        bitstream_length--;
        // printf("DEBUG: Loaded block %zu to main: %016lx\n", current_block++, main_reg);
    }
    if (bitstream_length > 0) {
        source_reg = *current_bitstream++;
        bitstream_length--;
        // printf("DEBUG: Loaded block %zu to source: %016lx\n", current_block++, source_reg);
    }

    size_t total_blocks = (sequence_length + 31) / 32;  // How many blocks we should have
    size_t total_number_of_blocks_left_to_move_to_source = total_blocks - 2;
    size_t remaining_bits_in_main = 64; // counts bits left in source and finally main reg

//asm("kmer_count_loop_start:");
    while (counted_kmers < total_kmers) {
        // printf("\nDEBUG: --- Iteration %zu ---\n", counted_kmers);
        // printf("main_reg: %016lx, source_reg: %016lx\n", main_reg, source_reg);
        // printf("remaining_bits_in_main: %zu\n", remaining_bits_in_main);

        // Extract and count kmer
        if (gap_length == 0) {
            kmer = main_reg & ((1ULL << (2 * kmer_length)) - 1);
        } else {
            uint64_t mask = (1ULL << (2 * (kmer_length - gap_position))) - 1;
            kmer = ((main_reg & mask) ^
                   ((main_reg >> (2 * gap_length)) & ~mask)) &
                   ((1ULL << (2 * kmer_length)) - 1);
        }

        // decode_kmer(kmer, kmer_length, gap_position, gap_length, kmer_str);
        // printf("DEBUG: Kmer %zu = %s from block %zu\n", counted_kmers, kmer_str, current_block - 1);
                             
        results[kmer]++;
        results[1ULL << (2 * kmer_length)]++;  // Increment total at the end of the array
        counted_kmers++;

        // Shift bits
        uint64_t new_bits = (source_reg & 0x3) << (64 - shift_amount);
        // printf("DEBUG: new_bits: %016lx\n", new_bits);

        #ifdef __x86_64__
        asm volatile (
            "shrdq $2, %1, %0"
            : "+r" (main_reg)
            : "r" (source_reg)
            : "cc"
        );
        #else
        main_reg = (main_reg >> shift_amount) | new_bits;
        #endif

        source_reg >>= shift_amount;
        remaining_bits_in_main -= 2;

        // printf("DEBUG: After shift - main_reg: %016lx, source_reg: %016lx, blocks remaining:%li\n", main_reg, source_reg, total_number_of_blocks_remaining);

        if (last_source_block_moved == 0 && remaining_bits_in_main == 0 && total_number_of_blocks_left_to_move_to_source == 0) {
            remaining_bits_in_main = (sequence_length % 32) * 2;
            last_source_block_moved = 1;
            // printf("DEBUG: last block in source. Remaining blocks, bits %li, %li\n", total_number_of_blocks_left_to_move_to_source, remaining_bits_in_main);
        }
        
        // Critical part: loading next block
        if (last_source_block_moved == 0 && remaining_bits_in_main == 0 && total_number_of_blocks_left_to_move_to_source > 0) {
                source_reg = *current_bitstream++;
                total_number_of_blocks_left_to_move_to_source--;
                remaining_bits_in_main = 64;
                current_block++;
                // printf("DEBUG: new block moved to source. Remaining blocks %li, %li\n", total_number_of_blocks_left_to_move_to_source, remaining_bits_in_main);
            }

    }
//asm("kmer_count_loop_end:");
    
    // free(kmer_str);
    // printf("DEBUG: Counted %zu kmers\n", counted_kmers);
}

uint64_t***** allocate_5d_results(int shortest_kmer, int longest_kmer, int broader_gaps) {
    uint64_t***** results = malloc(NUM_FILES * sizeof(uint64_t****));
    
    for (int f = 0; f < NUM_FILES; f++) {
        results[f] = malloc((longest_kmer + 1) * sizeof(uint64_t***));
        
        for (int k = shortest_kmer; k <= longest_kmer; k++) {
            if (k <= 0) continue;
            
            uint64_t total_kmers = 1ULL << (2 * k);
            
            // Calculate center positions
            int centers[4] = {-1, -1, -1, -1};  // Store up to 4 centers
            int num_centers = 0;
            
            if (broader_gaps && k > 3) {
                if (k % 2 == 0) {
                    // For even length, allocate 3 center positions
                    centers[0] = (k / 2) - 1;  // Left of center
                    centers[1] = k / 2;        // Center
                    centers[2] = (k / 2) + 1;  // Right of center
                    num_centers = 3;
                } else {
                    // For odd length, allocate 4 center positions
                    centers[0] = ((k - 1) / 2) - 1;  // Left of centers
                    centers[1] = (k - 1) / 2;        // First center
                    centers[2] = ((k - 1) / 2) + 1;  // Second center
                    centers[3] = ((k - 1) / 2) + 2;  // Right of centers
                    num_centers = 4;
                }
            } else {
                // Original behavior
                if (k % 2 == 0) {
                    centers[0] = k / 2;
                    num_centers = 1;
                } else {
                    centers[0] = (k - 1) / 2;
                    centers[1] = ((k - 1) / 2) + 1;
                    num_centers = 2;
                }
            }
            
            // Allocate array for gap positions
            results[f][k] = malloc((k + 1) * sizeof(uint64_t**));
            for (int g = 0; g <= k; g++) {
                results[f][k][g] = NULL;  // Initialize all to NULL
            }
            
            // Allocate position 0 for ungapped kmers
            results[f][k][0] = malloc(sizeof(uint64_t*));
            results[f][k][0][0] = calloc(total_kmers + 1, sizeof(uint64_t));  // +1 for total at the end
            if (!results[f][k][0][0]) {
                fprintf(stderr, "Failed to allocate memory for ungapped kmers at k=%d (size=%llu)\n",
                        k, total_kmers);
                exit(1);
            }
            
            // Allocate memory for each center position
            for (int i = 0; i < num_centers; i++) {
                int center = centers[i];
                if (center != -1) {
                    results[f][k][center] = malloc((MAX_GAP_LENGTH + 1) * sizeof(uint64_t*));
                    for (int l = 1; l <= MAX_GAP_LENGTH; l++) {
                        results[f][k][center][l] = calloc(total_kmers + 1, sizeof(uint64_t));  // +1 for total at the end
                        if (!results[f][k][center][l]) {
                            fprintf(stderr, "Failed to allocate memory for gapped kmers at k=%d, center=%d, gap=%d (size=%llu)\n",
                                    k, center, l, total_kmers);
                            exit(1);
                        }
                    }
                }
            }
        }
    }
    return results;
}

// Matching free function that only frees allocated memory
void free_5d_results(uint64_t***** results, int shortest_kmer, int longest_kmer) {
    if (!results) return;
    
    for (int f = 0; f < NUM_FILES; f++) {
        if (!results[f]) continue;
        
        for (int k = shortest_kmer; k <= longest_kmer; k++) {
            if (k <= 0 || !results[f][k]) continue;
            
            // Free ungapped kmer storage
            if (results[f][k][0]) {
                free(results[f][k][0][0]);
                free(results[f][k][0]);
            }
            
            // Calculate center positions
            int center1, center2;
            if (k % 2 == 0) {
                center1 = k / 2;
                center2 = -1;
            } else {
                center1 = (k - 1) / 2;
                center2 = center1 + 1;
            }
            
            // Free first center position
            if (results[f][k][center1]) {
                for (int l = 1; l <= MAX_GAP_LENGTH; l++) {
                    free(results[f][k][center1][l]);
                }
                free(results[f][k][center1]);
            }
            
            // Free second center position if it exists
            if (center2 != -1 && results[f][k][center2]) {
                for (int l = 1; l <= MAX_GAP_LENGTH; l++) {
                    free(results[f][k][center2][l]);
                }
                free(results[f][k][center2]);
            }
            
            free(results[f][k]);
        }
        free(results[f]);
    }
    free(results);
}

void init_buffer(LineBuffer* buffer, int shortest_kmer, int longest_kmer, int file_index, const char* filename) {
    buffer->file_index = file_index;
    buffer->finished = 0;
    buffer->write_pos = 0;
    buffer->read_pos = 0;
    buffer->count = 0;
    buffer->shortest_kmer = shortest_kmer;
    buffer->longest_kmer = longest_kmer;
    
    pthread_mutex_init(&buffer->mutex, NULL);
    pthread_cond_init(&buffer->not_full, NULL);
    pthread_cond_init(&buffer->not_empty, NULL);
    
    buffer->sequences = malloc(BUFFER_SIZE * sizeof(char*));
    buffer->seq_lengths = malloc(BUFFER_SIZE * sizeof(size_t));
    for (int i = 0; i < BUFFER_SIZE; i++) {
        buffer->sequences[i] = malloc(1024 * sizeof(char));  // Adjust size as needed
    }
    
    buffer->file = fopen(filename, "r");
    if (!buffer->file) {
        fprintf(stderr, "Failed to open file: %s\n", filename);
        exit(1);
    }
}

void destroy_buffer(LineBuffer* buffer) {
    for (int i = 0; i < BUFFER_SIZE; i++) {
        free(buffer->sequences[i]);
    }
    free(buffer->sequences);
    free(buffer->seq_lengths);
    
    pthread_mutex_destroy(&buffer->mutex);
    pthread_cond_destroy(&buffer->not_full);
    pthread_cond_destroy(&buffer->not_empty);
    
    fclose(buffer->file);
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

void print_kmer_statistics(uint64_t***** results, int file_number, int kmer_length,
                        int gap_position, int gap_length, uint64_t kmer) {
   // Original kmer count
    //return;
    
   uint64_t orig_count = results[file_number][kmer_length][gap_position][gap_length][kmer];
   
   // 1. Calculate shift ratios - only for ungapped kmers
   uint64_t left_count = 0;
   uint64_t right_count = 0;
   double shift_ratio = 0.0;

    if (kmer_length < 4) shift_ratio = -1;
    else if (orig_count > 0) {  // Handle both gapped and ungapped kmers
        uint64_t max_count = 0;
        
        // Print original kmer and count
        // printf("\nOriginal: ");
        // Kmerprint(kmer, kmer_length, gap_position, gap_length);
        // printf(" Count: %lu", orig_count);

        // --- LEFT SHIFT (drop rightmost base) ---
        uint64_t shifted_kmer_l = kmer >> 2;
        for (int base = 0; base < 4; base++) {
            uint64_t candidate =
                ((uint64_t)base << (2 * (kmer_length - 1)))  // new leftmost base
                | shifted_kmer_l;                            // the remaining (k-1) bases

            uint64_t count_candidate;
            if (gap_length == 0) {
                count_candidate = results[file_number][kmer_length][0][0][candidate];
            } else {
                // For left shift, gap position moves right by 1
                count_candidate = results[file_number][kmer_length][gap_position + 1][gap_length][candidate];
            }
            
            // Print left-shifted candidate
            // printf("\nLeft shift: ");
            // if (gap_length == 0) Kmerprint(candidate, kmer_length, 0, 0);
            // else Kmerprint(candidate, kmer_length, gap_position + 1, gap_length);
            // printf(" Count: %lu", count_candidate);
            
            if (count_candidate > max_count) {
                max_count = count_candidate;
            }
        }

        // --- RIGHT SHIFT (drop leftmost base) ---
        uint64_t shifted_kmer_r =
            (kmer & ((1ULL << (2ULL * (kmer_length - 1))) - 1ULL)) << 2;

        for (int base = 0; base < 4; base++) {
            uint64_t candidate = shifted_kmer_r | (uint64_t)base;
            
            uint64_t count_candidate;
            if (gap_length == 0) {
                count_candidate = results[file_number][kmer_length][0][0][candidate];
            } else {
                // For right shift, gap position moves left by 1
                count_candidate = results[file_number][kmer_length][gap_position - 1][gap_length][candidate];
            }
            
            // Print right-shifted candidate
            // printf("\nRight shift: ");
            // if (gap_length == 0) Kmerprint(candidate, kmer_length, 0, 0);
            // else Kmerprint(candidate, kmer_length, gap_position - 1, gap_length);
            // printf(" Count: %lu", count_candidate);
            
            if (count_candidate > max_count) {
                max_count = count_candidate;
            }
        }

        // printf("\nMax shifted count: %lu\n", max_count);
        shift_ratio = (double)max_count / (double)orig_count;
    }
   
   // 2. Calculate IC and position-specific variance
   double total_ic = 0.0;
   double total_var = 0.0;
   int valid_positions = 0;
   
   // For each position
   for(int pos = 0; pos < kmer_length; pos++) {
       // Get counts for consensus and mutations at this position
       uint64_t pos_counts[4] = {0};  // A,C,G,T counts
       
       int shift = 2 * pos;
       
       // Get counts for consensus and all mutations at this position
       uint64_t base_mask = ~(3ULL << shift);
       for(int base = 0; base < 4; base++) {
           uint64_t mut_kmer = (kmer & base_mask) | ((uint64_t)base << shift);
           pos_counts[base] = results[file_number][kmer_length][gap_position][gap_length][mut_kmer];
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
    if (shift_ratio < 0) {
        printf("\tNA\t%.3lf\t%.2lf%%\t%.2lf%%\t%.2lf%%",
               total_ic,
               100 * avg_variance/orig_count,
               100 * exp_variance/orig_count,
               100 * avg_variance/orig_count - 100 * exp_variance/orig_count);
    } else {
        printf("\t%.3lf\t%.3lf\t%.2lf%%\t%.2lf%%\t%.2lf%%",
               shift_ratio, total_ic,
               100 * avg_variance/orig_count,
               100 * exp_variance/orig_count,
               100 * avg_variance/orig_count - 100 * exp_variance/orig_count);
    }
}

// PRODUCER TO GENERATE LINE BUFFER FOR CONSUMERS
void* producer(void* arg) {
    LineBuffer* buffer = (LineBuffer*)arg;
    char line[1024];
    time_t last_print = time(NULL);
    
    while (fgets(line, sizeof(line), buffer->file)) {
        // Skip empty lines or lines with just whitespace
        size_t len = strlen(line);
        if (len == 0) continue;
        
        // Trim trailing whitespace and newlines
        while (len > 0 && (line[len-1] == '\n' || line[len-1] == '\r' || line[len-1] == ' ')) {
            line[--len] = '\0';
        }
        
        // Skip if line is empty after trimming
        if (len == 0) continue;
        
        pthread_mutex_lock(&buffer->mutex);
        
        while (buffer->count == BUFFER_SIZE) {
            time_t current = time(NULL);
            if (current - last_print >= 10) {
                // printf("Producer: Still waiting for buffer space after %ld seconds. Buffer count=%d\n",
                //       current - last_print, buffer->count);
                // fflush(stdout);
                last_print = current;
            }
            pthread_cond_wait(&buffer->not_full, &buffer->mutex);
        }
        
        strcpy(buffer->sequences[buffer->write_pos], line);
        buffer->seq_lengths[buffer->write_pos] = len;
        buffer->write_pos = (buffer->write_pos + 1) % BUFFER_SIZE;
        buffer->count++;
        
        pthread_cond_signal(&buffer->not_empty);
        pthread_mutex_unlock(&buffer->mutex);
    }
    
    pthread_mutex_lock(&buffer->mutex);
    buffer->finished = 1;
    pthread_cond_broadcast(&buffer->not_empty);
    pthread_mutex_unlock(&buffer->mutex);
    
    return NULL;
}

// CONSUMER FOR PWM GENERATION
// CONSUMER FOR PWM GENERATION
void* PWM_generator_consumer(void* arg) {
   PWMConsumerArgs* args = (PWMConsumerArgs*)arg;
   LineBuffer* buffer = args->buffer;
   size_t line_count = 1;
   struct match matches;
   struct count_pwm count_pwms[2];
   
   // initializes matches and count pwm for this file
   match_init(&matches, args->signal_pwm->width);
   count_pwm_init(&count_pwms[buffer->file_index], "pwm", args->signal_pwm->width, 0.0);

   // Calculate proper allocation size based on MAX_SEQ_LENGTH
   size_t num_blocks = (MAX_SEQ_LEN + 31) / 32;
   
   // Allocate thread-local bitstreams
   uint64_t* local_bitstream = (uint64_t*)calloc(num_blocks, sizeof(uint64_t));
   uint64_t* local_temp_bitstream = (uint64_t*)calloc(num_blocks, sizeof(uint64_t));
   if (!local_bitstream || !local_temp_bitstream) {
       fprintf(stderr, "Failed to allocate thread-local bitstreams\n");
       free(local_bitstream);
       free(local_temp_bitstream);
       return NULL;
   }

   // Get kmer length from args
   int kmer_length = args->kmer_length;

   while (1) {
       pthread_mutex_lock(&buffer->mutex);
       
       while (buffer->count == 0) {
           if (buffer->finished) {
               pthread_mutex_unlock(&buffer->mutex);
               
               // Calculate enrichment before returning - do this only once when processing is finished
               if (buffer->file_index == 1) { // Only do this once in signal file processing
                   size_t kmer_count_size = (1ULL << (2 * kmer_length));
                   
                   // Calculate total counts for normalization
                   uint64_t total_bg = 0;
                   uint64_t total_sig = 0;
                   
                   for (size_t i = 0; i < kmer_count_size; i++) {
                       total_bg += args->kmer_counts[0][i];
                       total_sig += args->kmer_counts[1][i];
                   }
                   
                   // Only proceed if we have counts
                   if (total_bg > 0 && total_sig > 0) {
                       // We need to sort arrays with the actual counts, not indices
                       uint64_t* bg_sorted = (uint64_t*)malloc(kmer_count_size * sizeof(uint64_t));
                       uint64_t* sig_sorted = (uint64_t*)malloc(kmer_count_size * sizeof(uint64_t));
                       
                       if (bg_sorted && sig_sorted) {
                           // Copy counts to sortable arrays
                           for (size_t i = 0; i < kmer_count_size; i++) {
                               bg_sorted[i] = args->kmer_counts[0][i];
                               sig_sorted[i] = args->kmer_counts[1][i];
                           }
                           
                           // Sort the arrays of counts directly
                           qsort(bg_sorted, kmer_count_size, sizeof(uint64_t), compare_counts);
                           qsort(sig_sorted, kmer_count_size, sizeof(uint64_t), compare_counts);
                           
                           // Calculate quartiles (each quartile is 25% of the data)
                           size_t quartile_size = kmer_count_size / 4;
                           
                           // Extract the two middle quartiles (25%-75%)
                           size_t middle_start = quartile_size;
                           size_t middle_end = 3 * quartile_size;
                           
                           // Sum the counts for the middle quartiles
                           uint64_t bg_middle_sum = 0;
                           uint64_t sig_middle_sum = 0;
                           for (size_t i = middle_start; i < middle_end; i++) {
                               bg_middle_sum += bg_sorted[i];
                               sig_middle_sum += sig_sorted[i];
                           }
                           
                           // Calculate lambda (protect against division by zero)
                           double lambda = 1.0; // Default to 1.0 (no correction) if we can't calculate
                           
                           if (args->background_output_pwm->enrichment) lambda = args->background_output_pwm->enrichment;
                           else {
                               if (bg_middle_sum > 0 && total_bg > 0 && total_sig > 0) {
                                   double bg_fraction = (double)bg_middle_sum / total_bg;
                                   double sig_fraction = (double)sig_middle_sum / total_sig;
                                   
                                   if (bg_fraction > 0) {
                                       lambda = sig_fraction / bg_fraction;
                                       //printf("\nLambda = %.2f, center fraction: %.2f, %.2f, middle %llu,%llu, total %llu,%llu", lambda, bg_fraction, sig_fraction, bg_middle_sum, sig_middle_sum, total_bg, total_sig);
                                   }
                               }
                           }
                           
                           // Simple safety check - if lambda is extreme, use a reasonable value
                           if (lambda <= 0.0 || lambda > 10.0) {
                               lambda = 1.0;
                           }
                           
                           // Calculate sizefactor from total kmer counts
                           // This matches the original code: sizefactor = number_of_sequences_analyzed[0] / number_of_sequences_analyzed[1]
                           // In our case, total_bg and total_sig are proportional to sequences
                           double sizefactor = 1.0; // Default
                           if (total_sig > 0) {
                               sizefactor = (double)total_bg / total_sig;
                               //lambda = (double) (sig_middle_sum * total_bg) / (double) (bg_middle_sum * total_sig);
                               //printf("\nLambda = %.2f", lambda);
                           }
                           
                           // Apply background correction to signal PWM exactly as in the provided code
                           if (args->signal_output_pwm && args->background_output_pwm) {
                               for (int i = 0; i < args->signal_output_pwm->width; i++) {
                                   for (int j = 0; j < 4; j++) {
                                       // Apply the exact formula from the provided code:
                                       // corrected = (signal * sizefactor) - (lambda * background)
                                       double signal_val = args->signal_output_pwm->incidence[j][i];
                                       double bg_val = args->background_output_pwm->incidence[j][i];
                                       
                                       double corrected = (signal_val * sizefactor) - (lambda * bg_val);
                                       // Zero out any negative values
                                       args->signal_output_pwm->incidence[j][i] = corrected > 0 ? corrected : 0;
                                   }
                               }
                           }
                           
                           // Store lambda in the PWM structure for reporting
                           if (args->signal_output_pwm) {
                               args->signal_output_pwm->enrichment=lambda;
                           }
                           if (args->background_output_pwm) {
                               args->background_output_pwm->enrichment=lambda;
                           }
                           
                           free(bg_sorted);
                           free(sig_sorted);
                       }
                   }
               }
               
               free(local_bitstream);
               free(local_temp_bitstream);
               return NULL;
           }
           pthread_cond_wait(&buffer->not_empty, &buffer->mutex);
           continue;  // Add this to recheck count after waking
       }
       
       char* sequence = buffer->sequences[buffer->read_pos];
       size_t seq_length = buffer->seq_lengths[buffer->read_pos];
       
       // Add validation HERE, before any processing
       if (!sequence || seq_length == 0) {
           pthread_mutex_unlock(&buffer->mutex);
           line_count++;
           continue;
       }
       
       buffer->read_pos = (buffer->read_pos + 1) % BUFFER_SIZE;
       buffer->count--;
       
       pthread_cond_signal(&buffer->not_full);
       
       size_t bitstream_length;
       memset(local_bitstream, 0, num_blocks * sizeof(uint64_t));
       
       if (encode_dna_sequence(sequence, local_bitstream, &bitstream_length) == 0) {
           pthread_mutex_unlock(&buffer->mutex);
           line_count++;
           continue;
       }
       
       pthread_mutex_unlock(&buffer->mutex);
       line_count++;

       __atomic_fetch_add((size_t*)&(args->total_reads[buffer->file_index]), 1, __ATOMIC_SEQ_CST);
       
       // Count kmers in forward strand for ALL reads (if sequence length is sufficient)
       if (seq_length >= kmer_length) {
           // Call count_kmers function with gap of 0 and gap length of 0
           count_kmers(args->kmer_counts[buffer->file_index],
                      local_bitstream,
                      seq_length,
                      kmer_length,
                      0,  // gap_position
                      0); // gap_length
       }
       
       // Create reverse complement for ALL reads
       reverse_complement_bitstream(local_bitstream, local_temp_bitstream, seq_length);
       
       // Count kmers in reverse complement strand for ALL reads
       if (seq_length >= kmer_length) {
           count_kmers(args->kmer_counts[buffer->file_index],
                      local_temp_bitstream,
                      seq_length,
                      kmer_length,
                      0,  // gap_position
                      0); // gap_length
       }
       
       // Check for PWM matches in forward strand
       if (Findpwmmatch_from_bitstream(args->signal_pwm, args->cutoff,
                                     local_bitstream, seq_length, &matches)) {
           
           // Increments match counts
           __atomic_fetch_add((size_t*)&(args->matched_reads[buffer->file_index]),
                             matches.position[0],
                             __ATOMIC_SEQ_CST);
           
           for (int i = 1; i <= matches.position[0]; i++) {
               struct count_pwm* output_pwm = (buffer->file_index == 0) ?
                   args->background_output_pwm : args->signal_output_pwm;
               
               Multinomial_add_to_pwm_from_bitstream(
                   output_pwm,
                   args->signal_pwm,
                   matches.position[i],
                   matches.score[i],
                   args->cutoff,
                   local_bitstream,
                   local_temp_bitstream,
                   seq_length,
                   args->background_pwm);
           }
       }
       
       // Check for PWM matches in reverse complement strand
       if (Findpwmmatch_from_bitstream(args->signal_pwm, args->cutoff,
                                     local_temp_bitstream, seq_length, &matches)) {
           
           __atomic_fetch_add((size_t*)&(args->matched_reads[buffer->file_index]),
                             matches.position[0],
                             __ATOMIC_SEQ_CST);
           
           for (int i = 1; i <= matches.position[0]; i++) {
               struct count_pwm* output_pwm = (buffer->file_index == 0) ?
                   args->background_output_pwm : args->signal_output_pwm;
               
               Multinomial_add_to_pwm_from_bitstream(
                   output_pwm,
                   args->signal_pwm,
                   matches.position[i],
                   matches.score[i],
                   args->cutoff,
                   local_temp_bitstream,
                   local_bitstream,  // Note: bitstreams are swapped here
                   seq_length,
                   args->background_pwm);
           }
       }
   }
   
   free(local_bitstream);
   free(local_temp_bitstream);
   return NULL;
}


// CONSUMER FOR KMER COUNTING
void* consumer(void* arg) {
    LineBuffer* buffer = (LineBuffer*)arg;
    size_t line_count = 1;
    
    while (1) {
        pthread_mutex_lock(&buffer->mutex);
        
        while (buffer->count == 0) {
            if (buffer->finished) {
                pthread_mutex_unlock(&buffer->mutex);
                return NULL;
            }
            pthread_cond_wait(&buffer->not_empty, &buffer->mutex);
        }
        
        char* sequence = buffer->sequences[buffer->read_pos];
        size_t seq_length = buffer->seq_lengths[buffer->read_pos];
        
        buffer->read_pos = (buffer->read_pos + 1) % BUFFER_SIZE;
        buffer->count--;
        
        pthread_cond_signal(&buffer->not_full);

        size_t bitstream_length;  // Added declaration here
        memset(buffer->bitstream, 0, 64 * sizeof(uint64_t));
        
        if (encode_dna_sequence(sequence, buffer->bitstream, &bitstream_length) == 0) {
            printf("\nError in DNA sequence on line %li, sequence rejected", line_count);
            pthread_mutex_unlock(&buffer->mutex);
            line_count++;
            continue;
        }
        
        pthread_mutex_unlock(&buffer->mutex);
        line_count++;
        
        // COUNT KMERS
        for(int strand = 0; strand < 2; strand++) {
            // Process ungapped kmers within valid range
            for (int k = buffer->shortest_kmer; k <= buffer->longest_kmer; k++) {
                if (k <= 0) continue;
                if (buffer->results[buffer->file_index][k][0] &&
                    buffer->results[buffer->file_index][k][0][0]) {
                    count_kmers(buffer->results[buffer->file_index][k][0][0],
                                buffer->bitstream, seq_length, k, 0, 0);
                }
            }
            
            // Process kmers within valid gap position and gap length range
            for (int k = buffer->shortest_kmer; k <= buffer->longest_kmer; k++) {
                if (k <= 0) continue;
                
                // Calculate gap positions based on kmer length and broader_gaps setting
                int gap_positions[4] = {-1, -1, -1, -1};  // Up to 4 positions
                int number_of_gap_positions;
                
                if (SHIFTED_GAP_POSITIONS && k > 3) {
                    if (k % 2 == 0) {
                        // For even length, use 3 positions
                        number_of_gap_positions = 3;
                        gap_positions[0] = (k/2) - 1;  // Left of center
                        gap_positions[1] = k/2;        // Center
                        gap_positions[2] = (k/2) + 1;  // Right of center
                    } else {
                        // For odd length, use 4 positions
                        number_of_gap_positions = 4;
                        gap_positions[0] = ((k-1)/2) - 1;  // Left of centers
                        gap_positions[1] = (k-1)/2;        // First center
                        gap_positions[2] = ((k-1)/2) + 1;  // Second center
                        gap_positions[3] = ((k-1)/2) + 2;  // Right of centers
                    }
                } else {
                    // Original behavior
                    if (k % 2 == 0) {
                        number_of_gap_positions = 1;
                        gap_positions[0] = k/2;
                    } else {
                        number_of_gap_positions = 2;
                        gap_positions[0] = k/2;
                        gap_positions[1] = k/2 + 1;
                    }
                }
                
                // For each gap position
                for (int i = 0; i < number_of_gap_positions; i++) {
                    int gap_pos = gap_positions[i];
                    if (gap_pos == -1) continue;  // Skip invalid positions
                    
                    // For each gap length (1-10)
                    for (int gap_len = 1; gap_len <= MAX_GAP_LENGTH; gap_len++) {
                        // Skip if kmer + gap would be longer than sequence
                        if (k + gap_len > seq_length) continue;
                        
                        if (buffer->results[buffer->file_index][k][gap_pos] &&
                            buffer->results[buffer->file_index][k][gap_pos][gap_len]) {
                            count_kmers(buffer->results[buffer->file_index][k][gap_pos][gap_len],
                                      buffer->bitstream, seq_length, k, gap_pos, gap_len);
                        }
                    }
                }
            }
            reverse_complement_bitstream(buffer->bitstream, buffer->temp_bitstream, seq_length);
        }
    }
}

// And modify print function to take file parameter:
void print_kmers_single_config(uint64_t***** results, int file, int kmer_len,
                             int shortest_kmer, int longest_kmer,
                             int gap_pos, int gap_len, double threshold,
                             int only_local_max, int background_correction, double length_diff_cutoff) {
    char* kmer_str = malloc(kmer_len + MAX_GAP_LENGTH + 1);
    if (!kmer_str) return;
    

    // if (gap_pos == 0 && gap_len == 0) {printf("UNGAPPED_CHECK k=%d first_kmer_count=%" PRIu64 "\n", kmer_len, results[file][kmer_len][0][0][0]);}
    
    uint64_t total_kmers = 1ULL << (2 * kmer_len);
    
    // printf("\nResults for k=%d gap_pos=%d gap_len=%d:\n",
    //       kmer_len, gap_pos, gap_len);
    
    for (uint64_t kmer = 0; kmer < total_kmers; kmer++) {
        uint64_t signal_count = results[1][kmer_len][gap_pos][gap_len][kmer];
        uint64_t background_count = results[0][kmer_len][gap_pos][gap_len][kmer];

        if (signal_count >= threshold) {
            short int is_localmax = Localmax((long int*****)results,
                                        1, kmer_len,
                                        shortest_kmer, longest_kmer,
                                        gap_pos, gap_len, kmer,
                                        length_diff_cutoff,
                                        background_correction,
                                        (double)results[1][kmer_len][gap_pos][gap_len][1ULL << (2 * kmer_len)] /
                                        (double)results[0][kmer_len][gap_pos][gap_len][1ULL << (2 * kmer_len)]);
            decode_kmer(kmer, kmer_len, gap_pos, gap_len, kmer_str);
            if (!only_local_max || is_localmax == 1) {
               printf("%s\t%i\t%i\t%" PRIu64 "\t%" PRIu64, kmer_str, gap_pos, gap_len, background_count, signal_count);
               if (is_localmax) {
                   printf("\tLocalmax");
                   // Add statistics for local maxima
                   print_kmer_statistics(results, 1, kmer_len, gap_pos, gap_len, kmer);
               }
               //printf("\n");
            }
            if (!only_local_max || is_localmax == 1) printf("\n");
            }
        }
    free(kmer_str);
}

int main(int argc, char* argv[]) {
    
    clock_t start = clock();
    
    short int max_seed_length = 20;
    short int count_kmers = 0;
    short int background_correction = 1;  // Default is to use background correction
    int kmer_length = 8;                 // Default kmer length
    double cmd_lambda = 0;              // Default: use calculated lambda value
    uint64_t***** results;
    double threshold;
    double length_diff_cutoff;
    int shortest_kmer;
    int longest_kmer;
    int print_only_localmax;
    struct count_pwm *multinomial_motif;
    short int multinomial;
    double lambda;
    char seed[max_seed_length+1];
    
    struct normalized_pwm qp;
    normalized_pwm_init(&qp, "pwm_for_multinomial", max_seed_length, 0);
    struct count_pwm count_pwms[2];
    
    if (argc < 4 || argc > 10) {
        fprintf(stderr, "Usage for kmer counting and local maxima detection: \n\t%s <background_file> <signal_file> <shortest_kmer_length> <longest_kmer_length> [count threshold for printing] [length_diff_cutoff]\n", argv[0]);
        fprintf(stderr, "\tlength_diff_cutoff: Optional parameter for setting the length difference cutoff for local maxima detection (default: 0.35 = longer kmer must have at least 35%% of the counts of the shorter one to be local max)\n");
        fprintf(stderr, "\nUsage for generating a PWM motif using a seed: \n\t%s <background_file> <signal_file> <IUPAC seed> <multinomial> [lambda] [-kl=<kmer_length>] [-noback] [-lambda=<lambda_value>]\n", argv[0]);
        fprintf(stderr, "\tmultinomial sets how many mismatches are allowed in kmers that contribute to motif\n");
        fprintf(stderr, "\t-kl=<kmer_length>: Optional parameter to set the k-mer length (default: 8)\n");
        fprintf(stderr, "\t-noback: Optional parameter to disable background correction\n");
        fprintf(stderr, "\t-lambda=<lambda_value>: Optional parameter to set lambda value manually\n");
        return 1;
    }

    // Parse optional arguments
    for (int i = 3; i < argc; i++) {
        if (strncmp(argv[i], "-kl=", 4) == 0) {
            kmer_length = atoi(argv[i] + 4);
        } else if (strcmp(argv[i], "-noback") == 0) {
            background_correction = 0;
        } else if (strncmp(argv[i], "-lambda=", 8) == 0) {
            cmd_lambda = atof(argv[i] + 8);
        }
    }

    char* filenames[NUM_FILES] = {argv[1], argv[2]};
    shortest_kmer = atoi(argv[3]);
    //printf("\nkmer length %i", base_kmer_length);
    
    // IF THIRD ARGUMENT IS A STRING, PARSES MOTIF GENERATION ARGUMENTS
    if (shortest_kmer == 0) {
        multinomial = atoi(argv[4]);
        lambda = (argc > 5 && argv[5][0] != '-') ? atof(argv[5]) : 0.0;
        strncpy(seed, argv[3], 20);
        Iupac_to_pwm(&qp, seed);
        //printf("\nMotif generation started with multinomial %i and lambda = %f", multinomial, lambda);
        //fflush(stdout);
        
        // Initialize the output PWMs
        for (int i = 0; i < 2; i++) {
            count_pwm_init(&count_pwms[i], qp.name, qp.width + 2 * FLANK_WIDTH + 1, 0.0);
        }
    }
    else {
        // ELSE PARSES LOCAL MAX ARGUMENTS AND COUNTS KMERS
        count_kmers = 1;
        longest_kmer = atoi(argv[4]);
        threshold = (argc > 5) ? atof(argv[5]) : 10;             // Use command line value if provided
        length_diff_cutoff = argc > 6 ? atof(argv[6]) : 0.35;                  // Use command line value if provided
        if (argc > 7) {
            if (argv[7][0] == 'c') background_correction = 1;
            if (argv[7][0] == 's') background_correction = 2;
        }
        print_only_localmax = 1;
        // Allocate 5D results array
        results = allocate_5d_results(shortest_kmer, longest_kmer, SHIFTED_GAP_POSITIONS);
        if (!results) {
            fprintf(stderr, "ERROR: Memory allocation failed for results array\n");
            return 1;
        }
    }
    
    // Print parameters
    // printf("DEBUG: Running with parameters:\n");
    // printf("DEBUG: Base kmer length: %d\n", base_kmer_length);
    // printf("DEBUG: Length difference cutoff: %f\n", length_diff_cutoff);
    // printf("DEBUG: Will analyze lengths %d, %d, and %d\n",
    //        base_kmer_length-1, base_kmer_length, base_kmer_length+1);
    
    GenerateMask ();
    PWMConsumerArgs pwm_args;
    
    // Process each file
    for (int file_idx = 0; file_idx < NUM_FILES; file_idx++) {
        LineBuffer buffer;
        buffer.results = results;
        buffer.bitstream = (uint64_t*)calloc(MAX_SEQ_LEN * 2, sizeof(uint64_t));
        buffer.temp_bitstream = (uint64_t*)calloc(MAX_SEQ_LEN * 2, sizeof(uint64_t));
        
        if (!buffer.bitstream || !buffer.temp_bitstream) {
            fprintf(stderr, "ERROR: Memory allocation failed for bitstream\n");
            free_5d_results(results, shortest_kmer, longest_kmer);  // Fixed call
            return 1;
        }

        init_buffer(&buffer, shortest_kmer, longest_kmer, file_idx, filenames[file_idx]);

        pthread_t producer_thread, consumer_thread;
        pthread_create(&producer_thread, NULL, producer, &buffer);
        if (count_kmers == 1) pthread_create(&consumer_thread, NULL, consumer, &buffer);
        else {
            // Create and initialize the PWMConsumerArgs structure
            
            pwm_args.buffer = &buffer;
            pwm_args.signal_pwm = &qp;
            pwm_args.background_pwm = NULL;  // Or appropriate background PWM if you have one
            pwm_args.signal_output_pwm = &count_pwms[1];
            pwm_args.background_output_pwm = &count_pwms[0];
            pwm_args.cutoff = 0 - (double) multinomial - 0.0001;

            // Initialize arrays for file counts
            pwm_args.num_files = 2;  // For background and signal
            if (buffer.file_index == 0) {
                pwm_args.total_reads = calloc(pwm_args.num_files, sizeof(_Atomic(size_t)));
                pwm_args.matched_reads = calloc(pwm_args.num_files, sizeof(_Atomic(size_t)));
                
                // Initialize k-mer counting parameters
                pwm_args.kmer_length = kmer_length;  // Use specified kmer length
                pwm_args.signal_pwm->enrichment = cmd_lambda;    // Pass command-line lambda value
                
                // Calculate size needed for k-mer counts arrays (4^k possible k-mers)
                size_t kmer_count_size = 1ULL << (2 * pwm_args.kmer_length); // 4^k
                
                // Allocate and zero-initialize k-mer count arrays for both files
                pwm_args.kmer_counts = (uint64_t**)malloc(pwm_args.num_files * sizeof(uint64_t*));
                if (!pwm_args.kmer_counts) {
                    fprintf(stderr, "Failed to allocate memory for k-mer counting array pointers\n");
                    return 1;
                }
                
                for (int i = 0; i < pwm_args.num_files; i++) {
                    pwm_args.kmer_counts[i] = (uint64_t*)calloc(kmer_count_size+1, sizeof(uint64_t));
                    if (!pwm_args.kmer_counts[i]) {
                        fprintf(stderr, "Failed to allocate memory for k-mer counting array %d\n", i);
                        // Clean up previously allocated memory
                        for (int j = 0; j < i; j++) {
                            free(pwm_args.kmer_counts[j]);
                        }
                        free(pwm_args.kmer_counts);
                        return 1;
                    }
                }
            }
            pthread_create(&consumer_thread, NULL, PWM_generator_consumer, &pwm_args);
        }

        pthread_join(producer_thread, NULL);
        pthread_join(consumer_thread, NULL);

        destroy_buffer(&buffer);
        free(buffer.bitstream);
        free(buffer.temp_bitstream);
    }
    
    if (count_kmers == 1) {
        int file = 0;
        // PRINT KMERS
        // for (int file = 0; file < NUM_FILES; file++) {
        //    printf("\nFile %d results:\n", file);
        printf("kmer\tgap position\tgap length\tbackground\tsignal\tLocalmax\tshift_ratio\tIC\tavg_var\texp_var\n");
        
        // First print ungapped kmers above threshold for all lengths
        for (int klen = shortest_kmer; klen <= longest_kmer; klen++) {
            if (klen <= 0) continue;
            print_kmers_single_config(results, file, klen, shortest_kmer, longest_kmer,
                                    0, 0, threshold, print_only_localmax, background_correction, length_diff_cutoff);
        }
        
        // Then print gapped local maxima
        for (int klen = shortest_kmer; klen <= longest_kmer; klen++) {
            if (klen <= 0) continue;
            
            // Calculate center positions
            if (klen % 2 == 0) {
                // Even length: one center position
                int gap_pos = klen/2;
                for (int gap_len = 1; gap_len <= 10; gap_len++) {
                    print_kmers_single_config(results, file, klen, shortest_kmer, longest_kmer,
                                            gap_pos, gap_len, threshold, print_only_localmax, background_correction, length_diff_cutoff);
                }
            } else {
                // Odd length: two center positions
                for (int gap_pos = klen/2; gap_pos <= klen/2 + 1; gap_pos++) {
                    for (int gap_len = 1; gap_len <= 10; gap_len++) {
                        print_kmers_single_config(results, file, klen, shortest_kmer, longest_kmer,
                                                gap_pos, gap_len, threshold, print_only_localmax, background_correction, length_diff_cutoff);
                    }
                }
            }
        }
    }
    else {
        int seed_length = strlen(seed);
        if (background_correction) {
            printf("\nMotif from all matches to seed in signal file (background corrected using lambda=%.4f)",
                  count_pwms[1].enrichment);
        } else {
            printf("\nMotif from all matches to seed in signal file (no background correction)");
        }
        //PrintMotif(&count_pwms[0], multinomial, seed_length);  // Print the background PWM
        PrintMotif(&count_pwms[1], multinomial, seed_length);  // Print the signal PWM
        printf("\nMatch statistics:");
        printf("\nBackground: %zu matched / %zu total (%.4f%%)",
              pwm_args.matched_reads[0], pwm_args.total_reads[0],
              100.0 * pwm_args.matched_reads[0] / pwm_args.total_reads[0]);
        printf("\nSignal: %zu matched / %zu total (%.4f%%)\n",
              pwm_args.matched_reads[1], pwm_args.total_reads[1],
              100.0 * pwm_args.matched_reads[1] / pwm_args.total_reads[1]);
    }
    
    // Cleanup
    if (count_kmers == 1) {
        free_5d_results(results, shortest_kmer, longest_kmer);
    } else {
        // Free k-mer counting resources
        if (pwm_args.kmer_counts) {
            for (int i = 0; i < pwm_args.num_files; i++) {
                if (pwm_args.kmer_counts[i]) {
                    free(pwm_args.kmer_counts[i]);
                }
            }
            free(pwm_args.kmer_counts);
        }
        free(pwm_args.total_reads);
        free(pwm_args.matched_reads);
    }
    
    clock_t end = clock();
    double cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("\nTotal execution time: %.2f seconds\n", cpu_time_used);
    
    return 0;
}
