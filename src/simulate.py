import msprime, tskit
from pathlib import Path
import subprocess, tempfile, io

# build a demography with K isolated populations have same ancestor
# model ID = OutOfAfricaExtendedNeandertalAdmixturePulse_3I21
# description: Three population out-of-Africa with an extended pulse of Neandertal admixture into Europeans
def year_to_gen(kya, gen_time=29):
    return int(round(kya * 1000.0 / gen_time))

def build_demography(ne):
    dem = msprime.Demography()

    # add populations
    dem.add_population(name="YRI", initial_size=ne)
    dem.add_population(name="CEU", initial_size=ne)
    dem.add_population(name="Neandertal", initial_size=ne)
    dem.add_population(name="Root", initial_size=ne)
    dem.add_population(name="Ancestral", initial_size=ne)

    # time point
    Tnsplit = year_to_gen(290) # split Neandertal and modern humans
    Tooa = year_to_gen(73.95)     # out of Africa
    Tnstart = year_to_gen(50)  # Neandertal admixture start
    Tnend = year_to_gen(30)    # Neandertal end admixture

    # population split events
    dem.add_population_split(
        time=Tooa,
        derived=["YRI", "CEU"],
        ancestral="Ancestral",
    )

    dem.add_population_split(
        time=Tnsplit,
        derived=['Neandertal', 'Ancestral'],
        ancestral="Root",
    )

    dem.add_symmetric_migration_rate_change(
        time=Tnstart,
        populations=["CEU", "Neandertal"],
        rate=0.029,
    )

    dem.add_symmetric_migration_rate_change(
        time=Tnend,
        populations=["CEU", "Neandertal"],
        rate=0.0,
    )

    print(f"Demography build")
    return dem


def simulate_vcf_bgzf(
    dem, N_perpop, mu, rec, seed, l, model,
    vcf_path, label_path
):
    # Set up the simulation parameters
    samples = [msprime.SampleSet(N_perpop, population=p, ploidy=2) for p in ["YRI", "CEU", "Neandertal"]]
    print(f"Samples → \n {samples}")

    dem.sort_events()
    #simulate
    ts_anc = msprime.sim_ancestry(
        samples=samples,
        demography=dem,
        sequence_length=l,
        recombination_rate=rec,
        model=model,
        random_seed=seed,
    )

    #mutation
    ts = msprime.sim_mutations(
        ts_anc, rate=mu, model=msprime.JC69(), random_seed=seed)
    
    #write bgzf vcf
    out_gz = Path(f"{vcf_path}.gz")  

    with open(out_gz, "wb") as fout:
        proc = subprocess.Popen(["bgzip", "-c"], stdin=subprocess.PIPE, stdout=fout)
        with io.TextIOWrapper(proc.stdin, encoding="utf-8") as pipe:
            ts.write_vcf(pipe)
        ret = proc.wait()
        if ret != 0:
            raise RuntimeError("bgzip failed while compressing VCF")
        
    # output labels
    with open(label_path, "w") as w:
        w.write("sample\tpopulation\n")
        for i, ind in enumerate(ts.individuals()):
            pop_id = ts.node(ind.nodes[0]).population   # 整数 id
            # 从 TreeSequence 里取这个种群的名字；如果没有名字就退回到 pop{id}
            pop_name = ts.population(pop_id).metadata.get("name", f"pop{pop_id}")
            w.write(f"tsk_{i}\t{pop_name}\n")


    print('simulation done.')
    