/* extract_fasta_multi_exon:
 *  
 * Takes a fasta file and list of intervals (one-per-line)
 * to extract in format:
 *   <defline_token> <n> <start_1>  <end_1> <start_2> <end_2> ... <start_n> <end_n>  +/- <optional-token>
 * OR
 *   <defline_token>
 * 
 * In first case, the <n> subsequences (referred to as 'pieces' in the
 * code) defined by <start_i..end_i> will be concatenated together in
 * the same output sequence.  This was implemented to allow multiple
 * exon CDS sequences to be created as a single sequence record.
 * For single 'exon' sequences simply specify <n> as '1' and include
 * the start and end positions as <start_1> and <end_1>.
 * For all lines which specify boundaries, the boundaries must be
 * valid within the sequence within [1..L] where L is the length
 * of the full sequence named <defline_token>. 
 * 
 * In the second case of line formatting, the full sequence called <defline_token>
 * will be output, on the + strand.
 *
 * If a fasta file name is not given, then the file is assumed to be <stdin>
 *
 * Optionally a single token can exist after the '+/-'. If one exists it will
 * be appended to the end of the new name for the extracted sequence.
 * 
 * Original author: Richa Argawal
 * Modified by:     Eric Nawrocki (modified to allow multiple pieces per interval)
 * 
 * To compile:
 *  'Optimized' for speed (not sure if 'optimization' is significant):
 *    gcc -O3 -o extract_fasta_multi_exon extract_fasta_multi_exon.c
 * 
 *  For debugging:
 *    gcc -Wall -g -o extract_fasta_multi_exon extract_fasta_multi_exon.c
 *
 * [Tue Mar  3 16:05:41 2015]
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* TRUE and FALSE may be already defined, */
/* let's not get in the way of other definitions. */
 
#if     defined(TRUE)
#undef TRUE 
#endif
 
#if     defined(FALSE)
#undef FALSE
#endif
 
typedef enum {FALSE = 0, TRUE = 1} bool; /* boolean type */

#define PRINT_CUTOFF 80 

#define MAXLINELEN 500000
#define MAXPIECES 150 /* max number of allowed exons */
#define UNKNOWN -2
#define WHITE_SPACE     " \n\r\t\v" 

#define PLUS   1
#define MINUS  2

typedef struct {
  char *name;               /* sequence name */
  int   npieces;            /* number of 'pieces' subsequences we will join in final output of the sequence */
  int   start;              /* start position of full region (same as pstart[0], first piece) */
  int   end;                /* end position of full region (same as pend[npieces-1], final piece) */
  /* note: for pstart and pend we hardcode the size (MAXPIECES) so we can easily use qsort() */
  int   pstartA[MAXPIECES]; /* [0..p..npieces-1], start position of each piece (e.g. exon) */
  int   pendA[MAXPIECES];   /* [0..p..npieces-1], end position of each piece (e.g. exon) */
  int   strand;             /* PLUS or MINUS, applies to all pieces */
  char *opttok;             /* optional token to append to end of the name we assign this sequence, read from input */
} interval_record;

void read_intervals(char *);
int  compare_entries(const void *, const void *);

void process_fasta(char *);
int  find_index(char *);
void print_fasta(int, char *, char *);

/* EPN added functions */
void dump_intervals(interval_record *interval_data, int interval_count);     /* for debugging only */
void free_intervals(interval_record *interval_data, int interval_count);     /* for cleanup only */

interval_record *interval_data;
int interval_count;

int main(int argc, char *argv[])
{
  if ((argc != 2) && (argc != 3))
    {
      fprintf(stderr,"extract_fa <interval_list> [<fa_file>]\n");
      exit(EXIT_FAILURE);
    }

  read_intervals(argv[1]);

  if (interval_count > 0)
    {
      qsort(interval_data, interval_count, sizeof(interval_record), compare_entries);

      if (argc == 2)
        process_fasta("stdin");
      else
        process_fasta(argv[2]);
    }

  /* clean up */
  free_intervals(interval_data, interval_count);

  exit(EXIT_SUCCESS);
}

void read_intervals(filename)
     char *filename;
{
  char descriptions[MAXLINELEN+1];
  FILE *infile;
  char *token;
  int   p; /* counter over pieces of an interval */
  int  np; /* number of pieces in current interval */

  infile = fopen(filename, "r"); 
  if (infile == NULL)
    { 
      fprintf(stderr, "Cannot open %s\n", filename); 
      exit(EXIT_FAILURE);
    }

  interval_count = 0;
  while(fgets(descriptions, MAXLINELEN, infile) != NULL) 
    interval_count++;

  if (interval_count == 0)
    {
      fclose(infile);
      interval_data = NULL;
      return;
    }
  else
    rewind(infile);

  interval_data = (interval_record *) malloc(interval_count * sizeof(interval_record));
  if (interval_data == NULL)
    {
      fprintf(stderr,"No space for %d contigs\n", interval_count);
      exit(EXIT_FAILURE);
    }

  /* Read each line of the interval file and store its info in interval_data[] */
  /* Example lines:
   * 
   * JRHK01000002.1 1 961815 964505 -
   * JOQW01000002.1 1 101803 104475 +
   * KN275973.1 2 226623 226774 226854 229725 -
   * KN275973.1 3 1 50 70 100 130 2000 +
   *
   * First token is sequence name, second is number of 'pieces' <np>, then 
   * comes <np*2> tokens, the start/end position of each piece. These must
   * be in ascending order. Final token is '+' or '-' specifying positive
   * or negative strand.
   */
  interval_count = 0;
  while(fgets(descriptions, MAXLINELEN, infile) != NULL) 
    {
      token = strtok(descriptions, WHITE_SPACE);
      if (token == NULL)
        {
          fprintf(stderr,"No contig on line %d of %s\n", interval_count+1, filename);
          exit(EXIT_FAILURE);
        }
        
      interval_data[interval_count].name = (char *) strdup(token);
      interval_data[interval_count].opttok = NULL; /* initialize optional extra token to NULL */
      if (interval_data[interval_count].name == NULL)
        {
          fprintf(stderr,"Can't store %s\n", token);
          exit(EXIT_FAILURE);
        }

      token = strtok(NULL, WHITE_SPACE);
      if (token == NULL)
        {
          interval_data[interval_count].strand     = PLUS;
          interval_data[interval_count].pstartA[0] = UNKNOWN;
          interval_data[interval_count].pendA[0]   = UNKNOWN;
          interval_data[interval_count].npieces    = 1;
          interval_data[interval_count].start      = UNKNOWN;
          interval_data[interval_count].end        = UNKNOWN;
        }
      else
        {
          interval_data[interval_count].npieces = atoi(token);
          np = interval_data[interval_count].npieces; /* for convenience only */

          /* Some error checking: */
          /* make sure we have a positive number of pieces */
          if(np < 1) { 
            fprintf(stderr,"Less than one piece specified (interval %d), this is not allowed\n", interval_count+1);
            exit(EXIT_FAILURE);
          }            
          /* make sure we haven't exceeded our hardcoded maximum */
          if(np > MAXPIECES) { 
            fprintf(stderr,"Maximum number of pieces exceeded (interval %d) %d > %d\n", interval_count+1, np, MAXPIECES);
            exit(EXIT_FAILURE);
          }            

          /* read start and end positions for each piece */
          for(p = 0; p < np; p++) { 
            /* read the start position for this piece */
            token = strtok(NULL, WHITE_SPACE);
            if (token == NULL)
              {
                fprintf(stderr,"No interval start for piece %d (interval %d) in file %s \n", p+1, interval_count+1, filename);
                exit(EXIT_FAILURE);
              }
            interval_data[interval_count].pstartA[p] = atoi(token);
            /* check to make sure the start position is positive if we're in a multi-piece interval */
            if((np > 1) &&  /* multi-piece interval */
               (interval_data[interval_count].pstartA[p] < 1)) { /* with a start value < 1 */
              fprintf(stderr, "Start position in multi-piece interval < 1; this is not allowed (%s)\n", interval_data[interval_count].name);
              exit(EXIT_FAILURE);
            }

            /* read the end position for this piece */
            token = strtok(NULL, WHITE_SPACE);
            if (token == NULL)
              {
                fprintf(stderr,"No interval end for piece %d (interval %d) in file %s \n", p+1, interval_count+1, filename);
                exit(EXIT_FAILURE);
              }
            interval_data[interval_count].pendA[p] = atoi(token);
            /* check to make sure the end position is positive */
            if(interval_data[interval_count].pendA[p] < 1) { /* end value < 1 */
              fprintf(stderr, "End position < 1; this is not allowed (%s)\n", interval_data[interval_count].name);
              exit(EXIT_FAILURE);
            }

            /* make sure that this piece's start <= end */
            if(interval_data[interval_count].pstartA[p] > interval_data[interval_count].pendA[p]) { 
              fprintf(stderr,"Interval %d, piece %d, start > end (%d > %d)\n", interval_count+1, p+1, 
                      interval_data[interval_count].pstartA[p], 
                      interval_data[interval_count].pendA[p]);
              exit(EXIT_FAILURE);
            }
            /* make sure that this piece's start > previous piece's end */
            if((p > 0) && (interval_data[interval_count].pendA[p-1] >= interval_data[interval_count].pstartA[p])) { 
              fprintf(stderr,"Interval %d, piece %d (%d..%d) does not come after piece %d (%d..%d)\n", interval_count+1, 
                      p+1, interval_data[interval_count].pstartA[p],   interval_data[interval_count].pendA[p],
                      p,   interval_data[interval_count].pstartA[p-1], interval_data[interval_count].pendA[p-1]);
              exit(EXIT_FAILURE);
            }
          }
          /* done reading all pieces, set interval 'start' and 'end' */
          interval_data[interval_count].start = interval_data[interval_count].pstartA[0];  /* start of first piece */
          interval_data[interval_count].end   = interval_data[interval_count].pendA[np-1]; /* end of final piece */

          /* read the strand for this interval (it will apply to all pieces) */
          token = strtok(NULL, WHITE_SPACE);
          if (token == NULL)
            {
              fprintf(stderr,"No interval strand for interval %d of %s\n", interval_count+1, filename);
              exit(EXIT_FAILURE);
            }

          if (strcmp(token, "+") == 0)
            interval_data[interval_count].strand = PLUS;
          else
            interval_data[interval_count].strand = MINUS;

          /* we allow one additional token, but its optional */
          token = strtok(NULL, WHITE_SPACE);
          if (token != NULL)
            {
              interval_data[interval_count].opttok = (char *) strdup(token);
              /* should be final token */
              token = strtok(NULL, WHITE_SPACE);
              if (token != NULL) {
                fprintf(stderr,"Extra token for interval %d of %s\n", interval_count+1, filename);
                exit(EXIT_FAILURE);
              }
            }
        }
      interval_count++;
    }
  fclose(infile);
}

int compare_entries(const void *f1, const void *f2)
{    
  interval_record *first;
  interval_record *second;
  int value;

  first = (interval_record *) f1;
  second = (interval_record *) f2;

  value = strcmp(first->name, second->name);
  if (value != 0)
    return(value);

  if (first->start < second->start)
    return(-1);
  if (first->start > second->start)
    return(1);

  if (first->end < second->end)
    return(-1);
  if (first->end > second->end)
    return(1);

  if (first->strand < second->strand)
    return(-1);
  if (first->strand > second->strand)
    return(1);

  return(0);
}

int find_index(name)
     char *name;
{
  int index, left, right;
  int temp;

  left = 0;
  right = interval_count - 1;

  while (left <= right)
    {
      index = (int)((left+right)/2);
      temp = strcmp(name, interval_data[index].name);
      if (temp == 0)
        {
          index--;
          while ((index >= 0) &&
                 (strcmp(name, interval_data[index].name) == 0))
            index--;
          return(index+1);
        }
      if (temp < 0)
        right = index-1;
      if (temp > 0)
        left = index+1;
    }

  return(UNKNOWN);
}

void process_fasta(filename)
     char *filename;
{
  char descriptions[MAXLINELEN+1];
  char copy_line[MAXLINELEN+1];
  char defline[MAXLINELEN+1];
  char name[MAXLINELEN+1];
  char *fasta, *current;
  FILE *infile;
  char *token;
  int max, line, interval_index;
  int len, current_length;
  int CHUNK = 1000000;

  if (strcmp(filename,"stdin") == 0)
    infile = stdin;
  else
    {
      infile = fopen(filename, "r"); 
      if (infile == NULL)
        { 
          fprintf(stderr, "Cannot open %s\n", filename); 
          exit(0);
        }
    }

  line = 0;
  interval_index = UNKNOWN;
  fasta = NULL;
  max = 0;
  current_length = 0;
  while(fgets(descriptions, MAXLINELEN, infile) != NULL) 
    {
      line++;
      strcpy(copy_line, descriptions);

      token = strtok(descriptions, WHITE_SPACE);
      if (token != NULL)
        {
          strcpy(name, token);

          if (name[0] == '>')
            {
              /* at new sequence, print previous (sub)sequence if necessary */
              if (interval_index != UNKNOWN)
                {
                  print_fasta(interval_index, fasta, defline);
                  if (fasta != NULL)
                    fasta[0] = '\0';
                  current_length = 0;
                }

              interval_index = find_index(&(name[1]));
              strcpy(defline, copy_line);
            }
          else
            {
              if (interval_index != UNKNOWN)
                {
                  copy_line[strlen(copy_line)-1] = '\0';

                  len = current_length + strlen(copy_line);
                  if (len < max) /* still have space in 'fasta' string */
                    strcpy(&(fasta[current_length]), copy_line);
                  else /* out of space 'fasta' string, realloc */
                    {
                      max = ((len/CHUNK)+1) * CHUNK;
                      current = (char *) malloc((max+1) * sizeof(char));
                      if (current == NULL)
                        {
                          fprintf(stderr,"No space\n");
                          exit(EXIT_FAILURE);
                        }

                      if (current_length > 0)
                        strcpy(current, fasta);
                      strcpy(&(current[current_length]), copy_line);
                      if(fasta != NULL) free(fasta);
                      fasta = current;
                    }
                  current_length = len;
                }
            }
        }
    }
  /* finished final sequence, print previous (sub)sequence if necessary */
  if (interval_index != UNKNOWN)
    print_fasta(interval_index, fasta, defline);

  if (fasta != NULL) free(fasta);

  if (strcmp(filename,"stdin") != 0)
    fclose(infile);
}

void print_fasta(given_index, given_fasta, given_defline)
     int given_index;
     char *given_fasta;
     char *given_defline;
{
  int index, strand, temp, count, interval_length;
  int np;     /* number of pieces for current interval */
  int p;      /* counter over pieces */
  int pp;     /* p prime, allows us to reverse traversal order over p for negative strand */
  int start;  /* start for current piece */
  int end;    /* end for current piece */

  if (given_fasta == NULL)
    {
      fprintf(stderr,"No sequence for %s\n", interval_data[given_index].name);
      exit(EXIT_FAILURE);
    }

  interval_length = strlen(given_fasta);
  for (index = given_index;
       ((index < interval_count) &&
        (strcmp(interval_data[index].name, interval_data[given_index].name) == 0));
       index++)
    {
      np     = interval_data[index].npieces;
      strand = interval_data[index].strand;

      if (interval_data[index].pstartA[0] == UNKNOWN) 
        {
          /* two sanity checks */
          if(np != 1) { 
            fprintf(stderr, "Problem parsing intervals, start set as unknown for multipiece interval (%s)\n", interval_data[index].name);
            exit(EXIT_FAILURE);
          }
          if(interval_data[index].pendA[0] != UNKNOWN) { 
            fprintf(stderr, "Problem parsing intervals, start set as unknown, but end is not: (%s)\n", interval_data[index].name);
            exit(EXIT_FAILURE);
          }

          interval_data[index].pstartA[0] = 1;
          interval_data[index].pendA[0]   = interval_length;
        }
      else if (interval_data[index].end > interval_length) { 
        fprintf(stderr, "End position exceeds sequence length (%d > %d) for sequence %s\n", interval_data[index].end, interval_length, interval_data[index].name);
        exit(EXIT_FAILURE);
      }        
      
      /* output defline */
      printf(">%s:", interval_data[index].name);
      for(p = 0; p < np; p++) { 
        printf("%s%d_%s%d:", 
               interval_data[index].pstartA[p] == 1 ? "<" : "", 
               interval_data[index].pstartA[p], 
               interval_data[index].pendA[p] == interval_length ? ">" : "", 
               interval_data[index].pendA[p]);
      }
      if (strand == PLUS) { 
        printf("+");
      }
      else {
        printf("-"); 
      }
      if(interval_data[given_index].opttok != NULL) { 
        printf(":%s\n", interval_data[index].opttok);
      }
      else { 
        printf("\n");
      }
      
      count = 0;
      for(p = 0; p < np; p++) { /* for each piece... */
        /* if we're on the opposite strand, we need to do this in the reverse order */
        pp = (strand == PLUS) ? p : (np - 1) - p; 
        
        /* Make input based on offset 0 */
        start = interval_data[index].pstartA[pp] - 1;
        end   = interval_data[index].pendA[pp] - 1;

        if (strand == PLUS)
          {
            for (temp = start; temp <= end; temp++)
              {
                count++;
                printf("%c", given_fasta[temp]);
                if ((count%PRINT_CUTOFF) == 0)
                  printf("\n");
              }
            if ((p == (np-1)) && (count%PRINT_CUTOFF) != 0)
              printf("\n");
          }
        else
          {
            for (temp = end ; temp >= start; temp--)
              {
                count++;
                switch(given_fasta[temp]) {
                case 'A': printf("T"); break;
                case 'C': printf("G"); break;
                case 'G': printf("C"); break;
                case 'T': printf("A"); break;
                case 'a': printf("t"); break;
                case 'c': printf("g"); break;
                case 'g': printf("c"); break;
                case 't': printf("a"); break;
                default: printf("N");
                }
                if ((count%PRINT_CUTOFF) == 0)
                  printf("\n");
              }
            if ((p == (np-1)) && (count%PRINT_CUTOFF) != 0)
              printf("\n");
          }
      } /* end of 'for(p = 0; p < np; p++)' */
    }
}

/* EPN added functions 
 * Tue Mar  3 09:59:25 2015
 */

void dump_intervals(interval_record *interval_data, int interval_count) { 
  int i, p;

  for(i = 0; i < interval_count; i++) { 
    printf("interval_data[%d]: %s %d pieces (%d..%d) strand: %d\n", i+1, interval_data[i].name, interval_data[i].npieces, interval_data[i].start, interval_data[i].end, interval_data[i].strand);
    for(p = 0; p < interval_data[i].npieces; p++) { 
      printf("\tpiece %d: %d..%d\n", p+1, interval_data[i].pstartA[p], interval_data[i].pendA[p]);
    }
    printf("\n");
  }
}

void free_intervals(interval_record *interval_data, int interval_count) { 
  int i;

  if(interval_data != NULL) { 
    for(i = 0; i < interval_count; i++) { 
      if(interval_data[i].name != NULL) { 
        free(interval_data[i].name); 
        interval_data[i].name = NULL;
      }
    }
  }
  free(interval_data);
  interval_data = NULL;

  return;
}
