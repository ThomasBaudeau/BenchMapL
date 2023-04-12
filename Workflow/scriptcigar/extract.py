from nucleotide import Nucleotide
import pysam 

def extract_reads(l, perfect):
    """
    Extract all the reads from several BAM files.

    Parameters
    ----------
    *args : str
        x BAM filename
    perfect : str
        Name of the BAM file containing the perfect alignment
        
    Return
    ------
    - 1/ Dictionnary storing all reads contained in the "perfect" BAM file.
        --> Keys : "Number from the name of the read"
        --> Values : Read corresponding to the number (pysam.AlignedSegment object)
    
    - 2/ Dictionnary storing all reads for each BAM file from the different mapping tools
        --> Keys : file1 , file2, ..., filen
        --> Values : Another dictionnary containing all reads from one file in the same format that the dictionnary create for the "perfect" file.
    
    - 3/ Length of the sequence reference (int)
    """
    reads_perfect = {}
    bamfile = pysam.AlignmentFile(perfect, "rb")
    for r in bamfile:
        reads_perfect[r.query_name] = r

    all_files = {}
    c = 1
    for arg in l:
        reads = {}
        bamfile = pysam.AlignmentFile(arg, "rb")
        for r in bamfile:
            if r.query_name.split("_")[2]!='human':
                reads['_'.join(r.query_name.split("_")[0:2])  ] = r
        all_files[f"file{c}"] = reads
        c += 1
    length= bamfile.lengths[0]
    return reads_perfect, all_files, length

def lnucl(length):
    """
    Create a list of n Nucleotide object representing the reference sequence (n being the length of the sequence).

    Parameter
    ---------
    length : int
        Length of the reference sequence

    Return
    ------
    list
        List of Nucleotide object of the same length that the reference sequence
    """
    l_nucl = [Nucleotide(pos) for pos in range(1, length+1)]
    return l_nucl

def coverage(reads, l_nucl, param):
    """
    Calculate the expected or observed coverage of the reads contained in one BAM file for all the position of the reference sequence
    - Expected coverage if the file put in parameter is the "perfect" file with the perfect alignment
    - Observed coverage if the file come from a mapping tool

    Parameters
    ----------
    reads : dict
        Dictionnary obtained by the extract_reads fonction containing all the reads for one file
            --> Keys : "Number from the name of the read"
            --> Values : Read corresponding to the number (pysam.AlignedSegment object)
    l_nucl : list
        List of n Nucleotide object representing the reference sequence (n being the length of the reference sequence)
    param : str
        Indicate which coverage to calculate
    """
    for i in range(1, len(l_nucl)+1): 
        for r in reads.values():  
            start = r.pos+1
            stop = start + len(r.seq) -1
            if i >= start and i <= stop: 
                    l_nucl[i-1].increase_count(1, param) 


def score(perfect, l_reads, l_nucl, length):
    """
    Calculate a score representing the quality of the alignment of the same read realised by different mapping tools
    compared to the "perfect" alignment of this read for all the position of the reference sequence.
    More a read will be correctly aligned at a position considered by the different mapping tools more the score for this position 
    will be near of 1, while more a read will be wrongly aligned, more the score will be near of 0.

    Parameters
    ----------
    perfect : pysam.AlignedSegment object
        Read n from the "perfect" dictionnary storing the reads coming from the "perfect" file.
    l_reads : list[pysam.AlignedSegment object]
        List of all the read n coming from the different mapping tools BAM file.
    l_nucl : list
        List of n Nucleotide object representing the reference sequence (n being the length of the reference sequence)
    """
    lperfect=perfect.get_reference_positions(full_length=True)
    for read in l_reads:
        if read.is_mapped and not read.is_supplementary and not read.is_secondary:
            lread=read.get_reference_positions(full_length=True)
            if len(lread)==len(lperfect):
                for idx,pos in enumerate(lperfect):
                    if pos == lread[idx]:
                        if pos:
                            l_nucl[pos].set_score(1)
                            l_nucl[pos].increase_count()
                    else:
                        if type(pos)==int:
                            l_nucl[pos].set_score(-1)
                            l_nucl[pos].increase_count()
                        else:
                            l_nucl[lread[idx]].set_score(-1)
            else:
                print('error in size')
        # if read.is_unmapped:
        #     for idx, pos in enumerate(lperfect):
        #         if pos:
        #             l_nucl[pos].set_score(0)

            


def scores(args, perfect=None):
    save = pysam.set_verbosity(0)
    if perfect:
        perfect, reads, length = extract_reads(args, perfect=perfect)
        l_nucl = lnucl(length)   

        for name in perfect.keys():
            read_n = []
            for f in reads.keys():
                read_n.append(reads[f][name])
            score(perfect[name], read_n, l_nucl, length)
        return l_nucl


