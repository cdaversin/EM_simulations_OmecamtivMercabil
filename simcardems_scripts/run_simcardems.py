# Run simulation with settings from config and custom stimulus domain

from pathlib import Path
import dolfin
import matplotlib.pyplot as plt
import simcardems
import numpy as np
import sys
import ufl_legacy
import typing


logger = simcardems.utils.getLogger(__name__)
here = Path(__file__).absolute().parent
outdir = here / "results"

# Final time
T=8500000

# Frequency at which with export (.xdmf) files for visualization
save_freq=10

# Cell model : "fully_coupled_ORdmm_Land", "fully_coupled_Tor_Land", "explicit_ORdmm_Land", "pureEP_ORdmm_Land".
cellmodel="fully_coupled_Tor_Land"

# Cell model initialization (.json file or None)
cell_init_file=None

# Disease state : "healthy", "hf"
disease_state="hf"

# PCL : Pacing cycle length (ms)
PCL=1000

# Restart from previous simulation : False / True
# if True : reload from a state file state.h5 stored in .<here>/results/state.h5
load_state=False

# Drug dosage from factors file (.json)
drug_factors_file = "drug_factors/OM/OM_Tor_0p0uM.json"

def stimulus_domain(mesh):
    marker = 1
    # Wave:
    subdomain = dolfin.CompiledSubDomain("x[0] < 1.0")
    domain = dolfin.MeshFunction("size_t", mesh, mesh.topology().dim())
    domain.set_all(0)
    subdomain.mark(domain, marker)

    return simcardems.geometry.StimulusDomain(domain=domain, marker=marker)

# New class for EMcoupling that can overwrite the built-in function
class EMCoupling(simcardems.models.fully_coupled_Tor_Land.EMCoupling):
    def __init__(
        self,
        geometry,
        **state_params,
    ) -> None:
        super().__init__(geometry=geometry, **state_params)

    def register_datacollector(self, collector) -> None:
        super().register_datacollector(collector=collector)

    def ep_to_coupling(self):
        super().ep_to_coupling()

def INaL(vs, parameters):
    (
        v,
        CaMKt,
        m,
        h,
        j,
        hp,
        jp,
        mL,
        hL,
        hLp,
        a,
        iF,
        iS,
        ap,
        iFp,
        iSp,
        d,
        ff,
        fs,
        fcaf,
        fcas,
        jca,
        ffp,
        fcafp,
        nca,
        nca_i,
        kr_c0,
        kr_c1,
        kr_c2,
        kr_o,
        kr_i,
        xs1,
        xs2,
        Jrelnp,
        Jrelp,
        ki,
        kss,
        cass,
        cansr,
        cajsr,
        XS,
        XW,
        CaTrpn,
        TmB,
        Cd,
        cai,
        nai,
        nass,
    ) = vs

    # Assign parameters
    scale_INaL = parameters["scale_INaL"]
    nao = parameters["nao"]
    F = parameters["F"]
    R = parameters["R"]
    T = parameters["T"]
    CaMKo = parameters["CaMKo"]
    KmCaM = parameters["KmCaM"]
    KmCaMK = parameters["KmCaMK"]

    # Drug factor
    scale_drug_INaL = parameters["scale_drug_INaL"]

    # Population factors
    scale_popu_GNaL = parameters["scale_popu_GNaL"]

    HF_scaling_CaMKa = parameters["HF_scaling_CaMKa"]
    HF_scaling_GNaL = parameters["HF_scaling_GNaL"]

    # Init return args

    # Expressions for the CaMKt component
    CaMKb = CaMKo * (1.0 - CaMKt) / (1.0 + KmCaM / cass)
    CaMKa = (CaMKb + CaMKt) * HF_scaling_CaMKa

    # Expressions for the reversal potentials component
    ENa = R * T * ufl.ln(nao / nai) / F

    # Expressions for the INaL component
    GNaL = 0.0279 * scale_INaL * scale_drug_INaL * scale_popu_GNaL * HF_scaling_GNaL
    fINaLp = 1.0 / (1.0 + KmCaMK / CaMKa)
    return (-ENa + v) * ((1.0 - fINaLp) * hL + fINaLp * hLp) * GNaL * mL
    
config = simcardems.Config(
    outdir=outdir,
    coupling_type=cellmodel,
    T=T,
    cell_init_file=cell_init_file,
    show_progress_bar=False,
    save_freq=save_freq,
    disease_state=disease_state,
    dt_mech=0.5,
    mech_threshold=0.05,
    PCL=PCL,
    geometry_path="geometries/slab.h5",
    geometry_schema_path="geometries/slab.json",
    load_state=load_state,
    drug_factors_file=drug_factors_file,
)

if config.load_state:
    logger.info("Load previously saved state")
    coupling = EMCoupling.from_state(
        path=here / "results/state.h5",
        drug_factors_file=config.drug_factors_file,
        popu_factors_file=config.popu_factors_file,
        disease_state=config.disease_state,
        PCL=config.PCL,
    )

else:
    logger.info("Create new state")
    geo = simcardems.geometry.load_geometry(
        mesh_path="geometries/slab.h5",
        stimulus_domain=stimulus_domain,
    )
    coupling = simcardems.models.em_model.setup_EM_model_from_config(
        config=config,
        geometry=geo,
    )

# Run SimCardEMS with the above configuration
runner = simcardems.Runner.from_models(config=config, coupling=coupling)
runner.solve(T=config.T, save_freq=config.save_freq, show_progress_bar=config.show_progress_bar)

