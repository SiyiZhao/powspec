#include "define.h"
#include "read_file.h"
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>
#include <sys/stat.h>

#define BUF     512     // Maximum number of charactors in filenames and lines.

static inline float parse_float_big(const uint8_t *src) {
  uint32_t x;
  float f;
  x = ((uint32_t) src[3] << 0) | ((uint32_t) src[2] << 8) |
      ((uint32_t) src[1] << 16) | ((uint32_t) src[0] << 24);
  memcpy(&f, &x, 4);
  return f;
}

static float parse_float_little(const uint8_t *src) {
  uint32_t x;
  float f;
  x = ((uint32_t) src[0] << 0) | ((uint32_t) src[1] << 8) |
      ((uint32_t) src[2] << 16) | ((uint32_t) src[3] << 24);
  memcpy(&f, &x, 4);
  return f;
}
    
static inline double parse_double_big(const uint8_t *src) {
  uint64_t x;
  double d;
  x = ((uint64_t) src[7] << 0) | ((uint64_t) src[6] << 8) |
      ((uint64_t) src[5] << 16) | ((uint64_t) src[4] << 24) |
      ((uint64_t) src[3] << 32) | ((uint64_t) src[2] << 40) |
      ((uint64_t) src[1] << 48) | ((uint64_t) src[0] << 56);
  memcpy(&d, &x, 8);
  return d;
}
  
static inline double parse_double_little(const uint8_t *src) {
  uint64_t x;
  double d;
  x = ((uint64_t) src[0] << 0) | ((uint64_t) src[1] << 8) |
      ((uint64_t) src[2] << 16) | ((uint64_t) src[3] << 24) |
      ((uint64_t) src[4] << 32) | ((uint64_t) src[5] << 40) |
      ((uint64_t) src[6] << 48) | ((uint64_t) src[7] << 56);
  memcpy(&d, &x, 8);
  return d;
}
  

/******************************************************************************
Function `read_bigfile`:
  Read a bigfile format catalog.
Arguments:
  * `fname`:    the filename of the input catalog;
  * `data`:     a pointer to the structure storing galaxy data;
  * `num`:      the number of galaxies;
  * `verb`:  0 for concise outputs, 1 for detailed outputs.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
int read_bigfile(const char *fname, DATA **data, size_t *num, const int verb) {
  FILE *fp;
  int i, k, Nf, byte;
  size_t *fsize, max, j, n;
  char path[BUF], line[BUF], value[3][BUF];
  char *bname, *end;
  uint8_t *chunk;

  if (verb) printf("\n  Header: `%s'.\n", fname);

  if (!(fp = fopen(fname, "r"))) {
    P_ERR("failed to open the header file.\n");
    return POWSPEC_ERR_FILE;
  }

  for (i = 0; i < 3; i++) {
    if (!(fgets(line, BUF, fp))) {
      P_ERR("failed to read line %d of the header file.\n", i + 1);
      return POWSPEC_ERR_FILE;
    }
    memset(value[i], 0, BUF);
    sscanf(line, "%s %s", path, value[i]);
  }

  sscanf(value[1], "%d", &i);
  if (i != 3) {
    P_ERR("the dimension of the data is not 3.\n");
    return POWSPEC_ERR_CATA;
  }
  sscanf(value[2], "%d", &Nf);
  if (Nf <= 0) {
    P_ERR("no data files found in the header.\n");
    return POWSPEC_ERR_CATA;
  }

  // Currently only support float and double numbers.
  if ((value[0][0] != '<' && value[0][0] != '>') || value[0][1] != 'f' ||
      (value[0][2] != '4' && value[0][2] != '8')) {
    P_ERR("unknown data format: %s", value[0]);
    return POWSPEC_ERR_CATA;
  }
  byte = value[0][2] - 48;

  MY_ALLOC(fsize, size_t, Nf, the number of records);
  MY_ALLOC(bname, char, Nf * BUF, the name of data files);
  memset(bname, 0, Nf * BUF);
  *num = max = 0;
  for (i = 0; i < Nf; i++) {
    if (!(fgets(line, BUF, fp))) {
      P_ERR("failed to read line %d of the header file.\n", i + 4);
      return POWSPEC_ERR_FILE;
    }
    if (sscanf(line, "%s %zu", bname + i * BUF, fsize + i) != 2) {
      P_ERR("failed to resolve line %d of the header file.\n", i + 4);
      return POWSPEC_ERR_FILE;
    }
    if ((end = strchr(bname + i * BUF, ':')) != NULL) *end = '\0';
    *num += fsize[i];
    if (max < fsize[i]) max = fsize[i];
  }
  fclose(fp);
  if (*num <= 0) {
    P_ERR("no data records found in the header.\n");
    return POWSPEC_ERR_CATA;
  }

  // Check the data files.
  strcpy(path, fname);
  if ((end = strrchr(path, '/')) != NULL) *(end + 1) = '\0';
  for (i = 0; i < Nf; i++) {
    if (bname[i * BUF] == '\0') {
      P_ERR("name not found for file %d.\n", i + 1);
      return POWSPEC_ERR_FILE;
    }
    sprintf(line, "%s%s", path, bname + i * BUF);
    if (access(line, R_OK)) {
      P_ERR("cannot read file `%s'.\n", line);
      return POWSPEC_ERR_FILE;
    }
  }

  if (verb) {
    printf("  %d files to be read with ", Nf);
    if (value[0][2] == '4') printf("single precision.\n");
    else if (value[0][2] == '8') printf("double precision.\n");
    printf("  Number of objects: %zu\n  Allocating memory ...", *num);
  }

  MY_ALLOC(chunk, uint8_t, max * 3 * byte, reading data files);
  MY_ALLOC(*data, DATA, *num, the tracers);
  if (verb) {
    printf("\r  ~ %.3g Mb memory allocated for the tracers.\n"
        "  Reading ...  0%%", sizeof(DATA) * (*num) / (1024.0 * 1024.0));
    fflush(stdout);
  }

  // Reading data.
  n = 0;
  for (i = 0; i < Nf; i++) {
    sprintf(line, "%s%s", path, bname + i * BUF);
    if (!(fp = fopen(line, "r"))) {
      P_ERR("failed to open file `%s'.\n", line);
      return POWSPEC_ERR_FILE;
    }

    if (fread(chunk, sizeof(int8_t) * fsize[i] * 3 * byte, 1, fp) != 1) {
      P_ERR("failed to read file `%s'.\n", line);
      return POWSPEC_ERR_FILE;
    }

    max = 0;
    if (value[0][0] == '>') {           // Big endian.
      if (byte == 4) {          // float32
        for (j = 0; j < fsize[i]; j++) {
          for (k = 0; k < 3; k++)
            (*data)[n].x[k] = (float) parse_float_big(chunk + max + k * 4);
          max += 12;
          n += 1;
        }
      }
      else {                    // float64
        for (j = 0; j < fsize[i]; j++) {
          for (k = 0; k < 3; k++)
            (*data)[n].x[k] = (double) parse_double_big(chunk + max + k * 8);
          max += 24;
          n += 1;
        }
      }
    }
    else {                              // Little endian.
      if (byte == 4) {          // float32
        for (j = 0; j < fsize[i]; j++) {
          for (k = 0; k < 3; k++)
            (*data)[n].x[k] = (float) parse_float_little(chunk + max + k * 4);
          max += 12;
          n += 1;
        }
      }
      else {                    // float64
        for (j = 0; j < fsize[i]; j++) {
          for (k = 0; k < 3; k++)
            (*data)[n].x[k] = (double) parse_double_little(chunk + max + k * 8);
          max += 24;
          n += 1;
        }
      }
    }

    fclose(fp);
    if (verb && i != Nf - 1) {
      printf("\b\b\b\b%3d%%", (i + 1) / Nf);
      fflush(stdout);
    }
  }

  free(fsize);
  free(bname);
  free(chunk);
  if (verb)
    printf("\r  %zu objects recorded from %d files.\n", n, Nf);

  return 0;
}


