# Bioinfotools exercise 3
# Put your code instead of the 'pass' statements
import random
from Bio import Align
from Bio.Seq import Seq
from Bio.Align import substitution_matrices
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO, AlignIO
#from Bio.Align.Applications import ClustalwCommandline
from collections import defaultdict


class matrix:
    def __init__(self, n_rows, n_cols):
        assert isinstance(n_rows, int)
        assert isinstance(n_cols, int)
        assert n_rows > 0
        assert n_cols > 0 
        self.mat = [[0 for i in range(n_cols)] for i in range(n_rows)]
        self.rows = n_rows
        self.cols = n_cols

    def get(self, i, j):
        return self.mat[i][j]

    def set(self, i, j, x):
        self.mat[i][j] = x


def upgma(D, seq_names_lst):
    joined = None
    while(len(seq_names_lst)>1):
        a,b = lowest_cell(D)
        joined = (seq_names_lst[a],seq_names_lst[b])
        del(seq_names_lst[b])
        seq_names_lst[a] = joined
        update_matrix(D, a, b)
    return joined

def lowest_cell(D):
    min = float("inf")
    x,y = -1, -1
    for i in range(D.rows):
        for j in range(D.cols):
            if D.get(i,j) < min and i != j:
                x = i
                y = j
                min  = D.get(i,j)
    if x>y:
        y,x = x,y
    return x,y

def update_matrix(D,a,b):
    for i in range(0, a):
        D.set(a,i,(D.get(a,i) + D.get(b,i))//2)
    for i in range(a+1,b):
        D.set(a,i,(D.get(a,i) + D.get(b,i))//2)
    for i in range(b+1,D.cols):
        D.set(a,i,(D.get(a,i) + D.get(b,i))//2)
    for j in range(D.cols):
        D.set(j,a,D.get(a,j))
    for k in range(D.rows):
        del(D.mat[k][b])
    del(D.mat[b])
    D.rows -= 1
    D.cols -= 1
    

def pairwise(seq1,seq2):
    aligner = Align.PairwiseAligner()
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -5
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    if len(seq1)==0 or len(seq2)==0:
        return 0
    alignments = aligner.align(seq1, seq2)
    return alignments.score


def globalpw_dist(seq_lst):
    return distance_matrix_calc(seq_lst,pairwise, True)
    
def distance_matrix_calc(seq_lst, calc_function, extra):
    # Calculate the pairwise similarity scores
    n = len(seq_lst)
    similarity_matrix = matrix(n,n)
     
    for i in range(n):
        for j in range(n):
            score = calc_function(seq_lst[i], seq_lst[j])
            similarity_matrix.set(i,j,score)
            similarity_matrix.set(j,i,score)
    if extra:
        # Determine S_max
        S_max = max(max(similarity_matrix.mat[i]) for i in range(similarity_matrix.rows))
        
        # Convert similarity scores to distances
        distance_matrix = matrix(n,n)
        for i in range(n):
            for j in range(n):
                distance_matrix.set(i,j,S_max-similarity_matrix.get(i,j ) + 1)
        return distance_matrix
    else:
        return similarity_matrix


def count_kmer(seq, k=3):
    assert isinstance(k,int)
    assert k > 0 
    kmer_d = {}
    for i in range(0,len(seq)-k+1):
        if len(seq) < i+k:
            break
        curr = seq[i:i+k]
        if curr.seq in kmer_d:
            kmer_d[curr.seq] += 1
        else:
            kmer_d[curr.seq] = 1
    return kmer_d

def kmer_distance_calculator(d_i,d_j):
    s = set(d_i.keys())
    for key in d_j.keys():
        s.add(key)
    score = 0
    for key in s:
        score += (d_i.get(key,0) - d_j.get(key,0))**2
    return score ** 0.5

def kmer_dist(seq_lst, k=3):
    dict_lst = [count_kmer(seq_lst[i],k) for i in range(len(seq_lst))]
    d_mat = distance_matrix_calc(dict_lst,kmer_distance_calculator, False)
    return d_mat

def read_fasta_file(path):
    seq_records = list(SeqIO.parse(path, "fasta"))
    return seq_records

def get_splits(tree,d,fst):
    if isinstance(tree, str) or (tree is None):
        return
    if fst:
        d[tree] = 0
    else: 
        d[tree] += 1
    if isinstance(tree[0], tuple):
        get_splits(tree[0], d,fst)
    if isinstance(tree[1], tuple):
        get_splits(tree[1], d,fst)
    
def T_sorter(tree):
    if isinstance(tree, str):
        return tree
    if isinstance(tree, tuple):
        if isinstance(tree[0], str) and isinstance(tree[1], str):
            return tuple(sorted([tree[0], tree[1]]))
        if isinstance(tree[0], str):
            return tuple([T_sorter(tree[1]),tree[0]])
        elif isinstance(tree[1], str):
            return tuple([T_sorter(tree[0]),tree[1]])
        else:
            return tuple([T_sorter(tree[0]), T_sorter(tree[1])])
        
        
        

def eval_dist(seq_lst, msa_aln_path, dist_func=globalpw_dist):
    msa_seq = read_fasta_file(msa_aln_path)
    D = dist_func(seq_lst)
    names = [seq.name for seq in seq_lst]
    T = upgma(D, names)
    T = T_sorter(T)
    od = {}
    get_splits(T, od, True)
    
    for _ in range(100):
        num_columns = len(msa_seq[0])
        selected_columns = random.sample(range(num_columns), num_columns // 2)
        new_sequences = []
        selected_columns.sort()
        for sequence in msa_seq:
            new_seq = "".join(sequence[i] for i in selected_columns)
            new_seq = new_seq.replace("-", "")
            new_sequences.append(SeqRecord(Seq(new_seq), id=sequence.id, description=sequence.description))
        
        n_names = [seq.id for seq in new_sequences]
        D_new = dist_func(new_sequences)
        T_new = upgma(D_new, n_names)
        T_new = T_sorter(T_new)
        
        nd = defaultdict(int)
        get_splits(T_new, nd,False)
        
        for key in od.keys():
            od[key] += nd[key]
    return od


    """
    :param seq_lst: list
        list of n sequences S1, S2, ..., Sn
    :param msa_aln_path: str
        ClustalW FASTA alignment file
    :param dist_func: a distance function name
    :return: dict
        keys are tuples representing each of the splits in the UPGMA tree T
        values are the counters of the splits in the random permutations
    """



if __name__ == '__main__':
    pass
