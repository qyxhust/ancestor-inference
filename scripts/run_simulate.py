import yaml
from pathlib import Path
import stdpopsim
import msprime
from src.simulate import simulate_vcf, write_vcf

def run_simulate(model = "normal"):
    cfg  = yaml.safe_load(Path("config/default.yaml").read_text())
    params = cfg['msprime']
    left = params['left']
    right = left + params['l']
    n_per_pop = params['N_perpop']
    ploidy = params['ploidy']
    seed = params['seed']
    msprime_model = params['model']
    model_name = {
        "normal": "OutOfAfricaExtendedNeandertalAdmixturePulse_3I21",
        "hard":   "AncientEurope_4A21",
        "admix":  "AmericanAdmixture_4B18",
    }


    species = stdpopsim.get_species("HomSap")
    model = species.get_demographic_model(model_name[model])

    contig = species.get_contig(
        "chr1",
        mutation_rate = model.mutation_rate,
        recombination_rate = model.recombination_rate,
        left=left,
        right=right,
    )
    
    engine = stdpopsim.get_engine("msprime")
    pop_ids = [pop.name for pop in model.populations]
    # construct a dictionary of {pop_id: n_per_pop}
    pop2n = {pid: n_per_pop for pid in pop_ids}
    # construct samples
    samples = model.get_sample_sets(pop2n, ploidy=ploidy)
    ts_anc = engine.simulate(
        demographic_model=model,
        contig=contig,
        samples=samples,
        msprime_model = msprime_model,
        seed=seed,
    )

    write_vcf(ts_anc, 
        vcf_path= cfg['project']['simulate_vcfs'],
        label_path= cfg['project']['simulate_label']
    )
    print("[OK] Simulation complete.")
