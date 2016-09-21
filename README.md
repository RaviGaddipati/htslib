HTSlib is an implementation of a unified C library for accessing common file
formats, such as [SAM, CRAM and VCF][1], used for high-throughput sequencing
data, and is the core library used by [samtools][2] and [bcftools][3].
HTSlib only depends on [zlib][4].
It is known to be compatible with gcc, g++ and clang.

HTSlib implements a generalized BAM index, with file extension `.csi`
(coordinate-sorted index). The HTSlib file reader first looks for the new index
and then for the old if the new index is absent.

This project also includes the popular tabix indexer, which indexes both `.tbi`
and `.csi` formats, and the bgzip compression utility.

[1]: http://samtools.github.io/hts-specs/
[2]: http://github.com/samtools/samtools
[3]: http://samtools.github.io/bcftools/
[4]: http://zlib.net/

### Building HTSlib

See [INSTALL](INSTALL) for complete details.
[Release tarballs][download] contain generated files that have not been
committed to this repository, so building the code from a Git repository
requires extra steps:

```sh
autoheader     # If using configure, generate the header template...
autoconf       # ...and configure script (or use autoreconf to do both)
./configure    # Optional, needed for choosing optional functionality
make
make install
```

[download]: http://www.htslib.org/download/

### Notes from Ravi

# FASTA Indexing

`#include "faidx.h"` 

**Building an index**
Given a FASTA filename `fn` , the index is built with
`int fai_build(const char *fn)` 
This creates the index `fn.fai` . To load the index, use:
`faidx_t *fai_load(const char *fn)` 
This opens `fn.fai` . 


### **The number of sequences:**

`int faidx_nseq(const faidx_t *fai)` 


### **The name of the i’th sequence:**

`const char *faidx_iseq(const faidx_t *fai, int i)` 


### **The sequence length:**

`int faidx_seq_len(const faidx_t *fai, const char *seq)` 


### **Get a subsequence:**

`char *fai_fetch(const faidx_t *fai, const char *reg, int *len)` 
where `reg`  is of the form `chr2:10,000–20,000` . `len`  indicates the length of the returned sequence, negative if error.
Alternate:
`char *faidx_fetch_seq(const faidx_t *fai, const char *c_name, int p_beg_i, int p_end_i, int *len)` 
Where the caller provides the parsed region.
The caller frees the returned subsequence.
 
The index is destroyed with
`void fai_destroy(faidx_t *fai)` 



## VCF/BCF Files

`#include "vcfutils.h"` 


### **File Setup**

VCF, BCF, SAM, BAM, CRAM, etc.. files are all opened with `htsFile *hts_open(char  *filename, char *mode)`  Where mode is `[rwa][bceguxz0–9]` :

- `r` Read
- `w` Write
- `a` Append

With `r`  the rest is automatically determined. Otherwise:

- `b` Binary
- `c` CRAM
- `g` gzip compressed
- `u` uncompressed
- `z` bgzf compressed
- `[0–9]` zlib compression level

 `bcf_open` , `vcf_open` , `bcf_close` , etc… are just macros for `hts_open` and `hts_close` .


### **VCF Header**

To load headers for VCF/BCF files, use
 `bcf_hdr_t *bcf_hdr_read(htsFile *fp)` 
  
Header functions:

- `const char **bcf_hdr_seqnames(const bcf_hdr_t *h, int *nseqs)` 
  - Get a list of sequence names (CHROM). Caller freed list, but not each name.
  - Returns an ID, get actual name with `bcf_hdr_id2name(hdr, ID)` 
- `bcf_hdr_nsamples(hdr)` Macro that gets the number of samples.
- `int bcf_hdr_set_samples(bcf_hdr_t *hdr, const char *samples, int is_file)` 
  - Only load a subset of samples, most efficient when reading BAM directly. Otherwise record is converted to BAM and then unpacked.
  - With `is_file = 0`,  `samples` is a CSV string of samples to load, by name.
  - With `is_file = 1` , include samples from a given filename.
  - `-`  includes all, `NULL`  excludes all.


### **Getting a new record**

Records are stored in a `bcf1_t` struct. Members are:

- `int rid` is CHROM col
- `int pos`  is POS col
- `int rlen`  is REF length
- `float qual`  QUAL col
- `int n_info, n_allele, n_fmt, n_sample` 
- `kstring_t shared, indiv`  **What is this?**
- `bcf_dec_t d` Decoded record, set with `bcf_unpack` 
### 

 `int bcf_read(htsFile *fp, const bcf_hdr_t *h, bcf1_t *v)` 
 The read then needs to be unpacked with 
 `int bcf_unpack(bcf1_t *b, int which)` 
 `which`  defines what is unpacked, where:

- `BCF_UN_STR == 0b1` unpacks CHROM, POS, ID, REF, ALT
- `BCF_UN_SHR == 0b111` is  all shared information
- `BCF_UN_FMT == 0b1000` Unpacks FORMAT and each sample.
- `BCF_UN_ALL == 0b1111` Unpacks everything



### **Sample FORMAT information**

`bcf_get_format_*(bcf_hdr_t *hdr, bcf1_t *line, char* tag,T *dst, T *ndst)` 
where `*, T` can be `int32, float, or char` and defines the type. `*` can also be `string` but is an abstracted, slower, version of `char` .

- `hdr` BCF header
- `line` bcf1_t record
- `tag` What values to get. `bcf_get_genotypes(hdr, line, dst, ndst)` is a shortcut for `tag == "GT"` 
- `dst`  is a pointer for results
- `ndst` number of results. For `GT` this will be twice as big, i.e. `dst[0,1]`  will belong to the first sample, `dst[2,3]` to the second, etc…

This function is restricted to `set_samples` , and allocates memory if needed. `unpack` is also called if it has not been called yet so it does not need to be explicitly called.


### **Sample INFO fields**

`bcf_get_info_int32(hdr,line,tag,dst,ndst)` 

- Returns negative value for some kind of tag error


### **Variant classification**

`int bcf_get_variant_types(bcf1_t *rec)` 

- Returns one of `VCF_REF == 0, VCF_SNP == 1, VCF_MNP == 2, VCF_INDEL == 4, VCF_OTHER == 8` 

Destroy the header with `bcf_hdr_destroy(bcf_hdr_t *h)` .


## SAM/BAM Files

`#include "sam.h"` 


### **File Setup**

Same as in VCF/BCF files.


### **Indexing**

To build an index:
`int sam_index_build(const char *fn, int min_shift)` 
Where `fn`  is the filename, and `min_shift`  is 0 to generate a BAI, positive to generate a CSI. Return snegative on error.

To load:
`hts_idx_t *sam_index_load(htsFile *fp, const char *fn)` 


### **BAM Header**
    typedef struct {
    int32_t n_targets; // Num ref seqs
    int32_t ignore_sam_err;
    uint32_t l_text; // Length of plain text header
    uint32_t *target_len; // Length of ref seqs
    int8_t *cigar_tab;
    char **target_name; // ref names
    char *text; // Header plain text
    void *sdict;
    } bam_hdr_t;

The header is read with
`bam_hdr_t *bam_hdr_read(BGZF *fp)` 
**Is BGZF same as htsfile?**


### **Reading a record**

A record is stored in:

    typedef struct {
    bam1_core_t core; // Alignment information
    int l_data, m_data; // Current and max length of data
    uint8_t *data; // All the variable length data, qname-cigar-seq-qual-aux
    #ifndef BAM_NO_ID
    uint64_t id;
    #endif
    } bam1_t;

where:

    typedef struct {
    int32_t tid; // chromosome ID
    int32_t pos; // 0 based left coord
    uint32_t bin:16, qual:8, l_qname:8; // What is this?
    uint32_t flag:16; // bitwise flags
    uint32_t n_cigar:16; // number of cigar ops
    int32_t l_qseq; // length of read
    int32_t mtid; // chrom ID of next read in template
    int32_t mpos; // 0 based coord of next read
    int32_t isize;
    } bam1_core_t;
    


### **Read a record**
    int bam_read1(
    BGZF *fp, // File pointer
    bam1_t *b // Resulting BAM record
    )

The endpos is given by:
`int32_t bam_endpos(const bam1_t *b)` 

To use an iterator through a region with an index,

    hts_itr_t *sam_itr_queryi(
    const hts_idx_t *idx, // Index pointer
    int tid, // Chrom name
    int beg, // Beg pos
    int end // End pos
    )

Then the results obtained with:

    int sam_itr_next(
    htsfp, // File pointer
    itr, // iterator
    r // record pointer
    )

The record is valid if the return value is greater than 0.

The header is destroyed with `void bam_hdr_destroy(bam_hdr_t *h)` 