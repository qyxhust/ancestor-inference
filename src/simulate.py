import msprime, tskit
from pathlib import Path
import subprocess, tempfile, io
import pysam

# build a demography with K isolated populations have same ancestor
# model 1: mornal
# model ID = OutOfAfricaExtendedNeandertalAdmixturePulse_3I21
#model 2: hard
# model ID = 
#model 3: admix
# model ID = AmericanAdmixture_4B18
# description: Three population out-of-Africa with an extended pulse of Neandertal admixture into Europeans
    
def write_vcf(ts, vcf_path, label_path):
    import numpy as np
    vcf_path = Path(vcf_path)
    vcf_path.parent.mkdir(parents=True, exist_ok=True)

    # ---------- Write bgzipped VCF ----------
    with open(vcf_path, "wb") as fout:
        proc = subprocess.Popen(["bgzip", "-c"], stdin=subprocess.PIPE, stdout=fout)
        with io.TextIOWrapper(proc.stdin, encoding="utf-8") as pipe:
            ts.write_vcf(pipe)
        ret = proc.wait()
        if ret != 0:
            raise RuntimeError("bgzip failed while compressing VCF")
        
    pysam.tabix_index(str(vcf_path), preset="vcf", force=True)


    # ---------- Compute ancestry proportions ----------
    pop_names = [p.metadata["name"] for p in ts.populations()]
    n_pops = len(pop_names)

    # 每个样本的 q 值： individuals × populations
    Q = np.zeros((ts.num_individuals, n_pops))

    # 遍历所有区间
    for tree in ts.trees():
        interval_len = tree.interval[1] - tree.interval[0]

        for ind_id, ind in enumerate(ts.individuals()):
            # diploid → 两个节点
            nodes = ind.nodes
            for node in nodes:
                anc = tree.roots if tree.parent(node) == tskit.NULL else tree.parent(node)

                # ancestry population = population of the MRCA
                pop_id = tree.population(node)   # 树上 node 所属的 population
                Q[ind_id, pop_id] += interval_len

    # 归一化（除以全基因组总长度）
    Q = Q / np.sum(Q, axis=1, keepdims=True)


    # ---------- Write label file ----------
    with open(label_path, "w") as w:
        # 写表头
        header = ["sample", "population"] + [f"q_{name}" for name in pop_names]
        w.write("\t".join(header) + "\n")

        for i, ind in enumerate(ts.individuals()):
            pop_id = ts.node(ind.nodes[0]).population
            pop_name = ts.population(pop_id).metadata["name"]

            row = [f"tsk_{i}", pop_name] + [f"{x:.5f}" for x in Q[i]]
            w.write("\t".join(row) + "\n")

    print("simulation done.")
