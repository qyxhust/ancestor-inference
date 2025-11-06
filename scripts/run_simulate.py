import yaml
from pathlib import Path
from src.simulate import build_demography, simulate_vcf_bgzf

def run_simulate():
    cfg  = yaml.safe_load(Path("config/default.yaml").read_text())

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
        vcf_path= cfg['project']['simulate_vcfs'],
        label_path= cfg['project']['simulate_label'],
    )
    print("[OK] Simulation complete.")
    
