
"""
@func read_paf(paf_fname)

Reads a PAF file of overlaps and returns a list of overlaps as
ordered tuples with the following structure:

    qn - query sequence name
    ql - query sequence length
    qb - query start coordinate
    qe - query end coordinate
    rc - 0 means query and target on the same strand, 1 means opposite
    tn - target sequence name
    tl - target sequence length
    tb - target start coordinate on the original strand
    te - target end coordinate on the original strand
"""

def read_paf(paf_fname):
    overlaps = []
    types = (str, int, int, int, lambda x: int(x[0] == '-'), str, int, int, int)
    for line in open(paf_fname, "r"):
        tokens = line.rstrip().split()[:9]
        overlaps.append(tuple(t(tok) for t, tok in zip(types, tokens)))
    return overlaps

"""
@func read_fasta(fasta_fname)

Reads a FASTA file of sequences and returns a list of
sequences and a list of names. The names and sequences
are ordered so that the same positions in each list correspond
to the same sequences in the FASTA.
"""
def read_fasta(fasta_fname):
    seqs, names, seq, name = [], [], [], ""
    for line in open(fasta_fname):
        if line[0] == ">":
            if len(seq) > 0:
                seqs.append("".join(seq))
                names.append(name)
            name = line.lstrip(">").split()[0].rstrip()
            seq = []
        else: seq.append(line.rstrip())
    if len(seq) > 0:
        seqs.append("".join(seq))
        names.append(name)
    return seqs, names

"""
Directed graphs are represented as dictionaries of dictionaries.
For example, if g is a an assembly graph, then g[u][v] = l means
that the weight of the edge (u -> v) is equal to l.
"""

"""
@func asg_delete_verts(asg, vdel)

Creates a new graph from `asg` except where all vertices
in set `vdel` are now isolated.
"""
def asg_delete_verts(asg, vdel):
    g = {u : {} for u in range(len(asg))}
    for u in asg:
        if u in vdel: continue
        for v in asg[u]:
            if v in vdel: continue
            g[u][v] = asg[u][v]
    return g

"""
@func asg_delete_arcs(asg, adel)

Creates a new graph from `asg` except where all edges
stored by `adel` are removed. A deleted edge (u -> v) in the
dictionary `adel` (whose keys are source vertices and whose
values are sets of vertices) is represented by the existence of
v in the set adel[u].
"""
def asg_delete_arcs(asg, adel):
    g = {u : {} for u in range(len(asg))}
    for u in asg:
        for v in asg[u]:
            if not v in adel[u]:
                g[u][v] = asg[u][v]
    return g

"""
@func asg_delete_asymm(asg)

Creates a new graph from `asg` except where all asymmetric
edges are removed. A an edge (u -> v) is considered asymmetric
if the edge (comp(v) -> comp(u)) is not also in the graph.
"""
def asg_delete_asymm(asg):
    g = {u : {} for u in range(len(asg))}
    for u in asg:
        for v in asg[u]:
            uc, vc = u^1, v^1
            if uc in asg[vc]:
                g[u][v] = asg[u][v]
    return g

"""
@func gen_asm_graph(overlaps, seqs, names, max_hang, max_ratio)

Creates an assembly graph in the vein of miniasm. Takes a list of
overlaps, read sequences, and read names as inputs.
"""
def gen_asm_graph(overlaps, seqs, names, max_hang=1000, max_ratio=0.8):
    n = len(seqs)
    namemap = {names[i] : i for i in range(n)}
    asg = {i : {} for i in range(2*n)}
    vdel = set()
    for overlap in overlaps:
        qn, ql, qb, qe, rc, tn, tl, tb, te = overlap
        if qn == tn: continue
        qi, ti = namemap[qn], namemap[tn]
        assert ql == len(seqs[qi]) and tl == len(seqs[ti])
        b1, e1, l1 = qb, qe, ql
        b2, e2, l2 = (tb, te, tl) if rc == 0 else (tl - te, tl - tb, tl)
        overhang = min(b1, b2) + min(l1 - e1, l2 - e2)
        maplen = max(e1 - b1, e2 - b2)
        qf, tf = 2*qi, 2*ti if rc == 0 else 2*ti + 1
        if overhang > min(max_hang, maplen * max_ratio):
            continue
        elif b1 <= b2 and l1 - e1 <= l2 - e2:
            vdel.add(qf)
            vdel.add(qf^1)
            continue
        elif b1 >= b2 and l1 - e1 >= l2 - e2:
            vdel.add(tf)
            vdel.add(tf^1)
            continue
        elif b1 > b2:
            asg[qf][tf] = b1 - b2
            asg[tf^1][qf^1] = (l2 - e2) - (l1 - e1)
        else:
            asg[tf][qf] = b2 - b1
            asg[qf^1][tf^1] = (l1 - e1) - (l2 - e2)

    asg = asg_delete_verts(asg, vdel)
    asg = asg_delete_asymm(asg)
    return asg

"""
@func asm_graph_del_trans(asg, fuzz)

Uses Meyers (2005) linear-time transitive reduction algorithm.
"""
def asm_graph_del_trans(asg, fuzz):
    n = len(asg)
    vdel = set() # deleted vertices
    adel = {v : set() for v in range(n)} # deleted edges
    mark = [0] * n # all vertices start marked 0
    n_reduced = 0
    for v in range(n): # for each vertex v
        nv = len(asg[v]) # nv is outdegree
        if nv == 0: continue
        if v in vdel: # if v is to be deleted, then delete all outgoing edges
            for w in asg[v]:
                adel[v].add(w)
                n_reduced += 1
            continue
        for w in asg[v]: # mark all vertices reached from v as 1
            mark[w] = 1
        av = sorted(asg[v].items(), key=lambda x: x[1]) # av is array of outgoing edges from v, sorted by prefix length
        L = av[-1][1] + fuzz # get largest outgoing edge length + fuzz
        for w, vwl in av: # for each outgoing arc (v -> w)
            nw = len(asg[w]) # nw is outdegree of w
            if mark[w] != 1: continue # if w hasn't been marked as 1, continue
            aw = sorted(asg[w].items(), key=lambda x: x[1]) # aw is array of outgoing edges from w, sorted by prefix length
            for x, wxl in aw: # for each outgoing arc w ->x
                if wxl + vwl > L: break # if sum of transitive prefixes is greater than L, then stop
                if mark[x] != 0: mark[x] = 2 # if x is not marked 0, then mark it as 2
        for w, _ in av: # for each outgoing arc (v -> w)
            if mark[w] == 2: # if v is marked 2, then delete (v->w)
                adel[v].add(w)
                n_reduced += 1
            mark[w] = 0 # mark w as 0.

    asg = asg_delete_verts(asg, vdel)
    asg = asg_delete_arcs(asg, adel)
    asg = asg_delete_asymm(asg)
    return asg
