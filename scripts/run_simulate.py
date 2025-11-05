import yaml
from pathlib import Path
from src.simulate import build_demography, simulate_vcf_bgzf

def run_simulate():
    cfg  = yaml.safe_load(Path("config/default.yaml").read_text())
    outdir = cfg['project']['outdir']

    params = cfg['msprime']
    dem = build_demography(
        ne=params['ne'],
    )
    simulate_vcf_bgzf(
        dem=dem,
        N_perpop=params['N_perpop'],
        mu=params['mu'],
        rec=params['rec'],
        seed=params['seed'],
        l=params['l'],
        model="hudson",
        vcf_path=Path(outdir) / "simulate" / "truth.vcf",
        label_path=Path(outdir) / "simulate" / "labels.tsv",
    )
    print("[OK] Simulation complete.")
    
