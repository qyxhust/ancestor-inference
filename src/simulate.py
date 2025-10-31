import msprime, tskit
from pathlib import Path
import subprocess, tempfile, io

# build a demography with K isolated populations have same ancestor
def build_demography(k, ne):
    dem = msprime.Demography()
    dem.add_population(name="ancestral", initial_size=ne)
    for p in range(k):
        dem.add_population(name=f"pop{p}", initial_size=ne)
    dem.add_population_split(
        time=5000,
        derived=[f"pop{i}" for i in range(k)],
        ancestral="ancestral",
    )
    print(f"Demography build")
    return dem

# 生成长度为 L、全部为 'A' 的参考 FASTA，并用 samtools 建索引
from pathlib import Path
import subprocess

def write_allA_ref(L, ref_path: Path, chrom="1"):
    """
    L: 参考序列长度（int，如 100_000）
    ref_path: 输出 FASTA 路径（Path，如 Path("data/simulate/ref.fa")）
    chrom: FASTA 头名（需与 VCF 的 CHROM 一致，默认 '1'）
    """
    ref_path.parent.mkdir(parents=True, exist_ok=True)

    # 按 60 列换行写入（常见 FASTA 风格）
    with open(ref_path, "w") as f:
        f.write(f">{chrom}\n")
        row = "A" * 60
        full_len = int(L)
        # 写满 L 个 'A'
        for i in range(0, full_len, 60):
            end = i + 60
            if end <= full_len:
                f.write(row + "\n")
            else:
                f.write("A" * (full_len - i) + "\n")


def simulate_vcf_bgzf(
    dem, k, N_perpop, mu, rec, seed, l, model,
    vcf_path, label_path, ref_path
):
    # Set up the simulation parameters
    samples = [msprime.SampleSet(N_perpop, population=p, ploidy=2) for p in range(k)]
    print(f"Samples → \n {samples}")

    #simulate
    ts_anc = msprime.sim_ancestry(
        samples=samples,
        demography=dem,
        sequence_length=l,
        recombination_rate=rec,
        model=model,
        random_seed=42,
    )

    #mutation
    ts = msprime.sim_mutations(
        ts_anc, rate=mu, model=msprime.JC69(), random_seed=seed)
    
    #write reference fasta
    write_allA_ref(l, ref_path, chrom="1")
    
    #write bgzf vcf
    out_gz = Path(f"{vcf_path}.gz")  

    with open(out_gz, "wb") as fout:
        proc = subprocess.Popen(["bgzip", "-c"], stdin=subprocess.PIPE, stdout=fout)
        with io.TextIOWrapper(proc.stdin, encoding="utf-8") as pipe:
            ts.write_vcf(pipe)
        ret = proc.wait()
        if ret != 0:
            raise RuntimeError("bgzip failed while compressing VCF")
        
    #output labels
    with open(label_path, "w") as w:
        w.write("sample\tpopulation\n")
        for i, ind in enumerate(ts.individuals()):
            pop = ts.node(ind.nodes[0]).population
            w.write(f"tsk_{i}\tpop{pop}\n")

    print('simulation done.')
    