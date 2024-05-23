#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  2 21:53:41 2023

@author: shaun
"""
import pyduo
from pyduo.containers import MoleculeInformation, Grid, DiagonalizerBlock
from pyduo.electronic import ElectronicState
from pyduo.molecular import MolecularSystem, OffDiagonalComponents
from pyduo.components import PotentialEnergy, BobRot, SpinRot, SpinOrbit, Lambda_q, Lambda_p2q, Lambda_opq, SpinOrbitX
from pyduo.functionals import EMO, BOBLEROY, POLYNOM_DECAY_24, GRID, LORENTZ, MLR, POLYNOM
from pyduo.api import process_blocks

import glob

from pyduo.api import pot_file_dataframes, plot_abinitio_calc

import pandas as pd
    

#%%
prev = process_blocks('./SrH_v372.inp')
prev_X2Sig = {'poten': prev['poten']['X2Sig+'],
              'bob-rot': prev['bob-rot']['<X2Sig+|BR|X2Sig+>'],
              'spin-rot': prev['spin-rot']['<X2Sig+|SR|X2Sig+>']}



    
#%% DUO Initialisation commands

molecule_info = MoleculeInformation("Sr", "H", 87.905612, 1.007825035)

grid_block = Grid(861, 1.40, 10.00, 0)

system = MolecularSystem()

system.molecule("Sr", "H", 87.905612, 1.007825035)
    
system.integration_grid(grid_block)

system.set_solution_method("SINC")

system.diagonalizer = DiagonalizerBlock()
system.initialise()

#%% X2Sigma+ Setup
electronic_state = ElectronicState('X2Sig+')


block = PotentialEnergy(1, 'X2Sig+', 0, 2, EMO(**prev_X2Sig['poten']),symmetry="+")
block.set_ab_initio_data(pd.read_csv('./PEC_X2SIG.csv'))
electronic_state.potential(block)

bobl_functional = BOBLEROY(**prev_X2Sig['bob-rot'])
bob_rot = BobRot(bobl_functional,potential_energy_instance=block)
electronic_state.bobrot(bob_rot)

spinrot_functional = BOBLEROY(**prev_X2Sig['spin-rot'])
spin_rot = SpinRot(spinrot_functional, potential_energy_instance=block)
electronic_state.spinrot(spin_rot)

system.add_state(electronic_state)

#%% E2Pi Setup
try:
    
    prev_E2Pi = {'poten': prev['poten']['E2Pi'],
                  'bob-rot': prev['bob-rot']['<E2Pi|BR|E2Pi>'],
                  'spin-rot': prev['spin-rot']['<E2Pi|SR|E2Pi>'],
                  'lambda-q': prev['lambda-q']['<E2Pi|lambda-q|E2Pi>'],
                  'lambda-p2q': prev['lambda-p2q']['<E2Pi|lambda-p2q|E2Pi>'],
                  'spin-orbit': prev['spin-orbit-x']['<E2Pi|LS|E2Pi>']}
    
    E2Pi_LSZ = {'RE': 2.1083727,
                 'P': 1,
                 'A_values': [0.00931721465064401,
                  0.729904639549337],
                 'AINF': 1}
    
    E2Pi = ElectronicState('E2Pi')
    
    E2Pi.potential(PotentialEnergy(2, "E2Pi", mult=2, lambda_=1, functional_instance = EMO(**prev_E2Pi['poten'])))
    E2Pi.potential_energy.set_ab_initio_data(pd.read_csv('./SrH_E_Poten.csv'))
    
    E2Pi.spinorbitcoupling(SpinOrbitX(functional_instance = POLYNOM_DECAY_24(**prev_E2Pi['spin-orbit']), potential_energy_instance=E2Pi.potential_energy))
    E2Pi.spin_orbit_coupling.set_ab_initio_data(pd.read_csv('./LSZ_E_E.csv'))
    
    E2Pi.spinrot(SpinRot(BOBLEROY(**prev_E2Pi['spin-rot']),potential_energy_instance=E2Pi.potential_energy))
    E2Pi.bobrot(BobRot(BOBLEROY(**prev_E2Pi['bob-rot']),potential_energy_instance=E2Pi.potential_energy))
    
    E2Pi.lambdap2q(Lambda_p2q(BOBLEROY(**prev_E2Pi['lambda-p2q']), potential_energy_instance=E2Pi.potential_energy))
    E2Pi.lambdaq(Lambda_q(BOBLEROY(**prev_E2Pi['lambda-q']), potential_energy_instance=E2Pi.potential_energy))
    
    system.add_state(E2Pi)
except:
    pass
#%%

prev_A2Pi = {'poten': prev['poten']['A2Pi'],
              'bob-rot': prev['bob-rot']['<A2Pi|BR|A2Pi>'],
              'spin-rot': prev['spin-rot']['<A2Pi|SR|A2Pi>'],
              'lambda-q': prev['lambda-q']['<A2Pi|lambda-q|A2Pi>'],
              'lambda-p2q': prev['lambda-p2q']['<A2Pi|lambda-p2q|A2Pi>'],
              'spin-orbit': prev['spin-orbit-x']['<A2Pi|LS|A2Pi>']}

A2Pi_Potential_dict = {"V0": 13500.571,
                       "RE": 2.1209,#86Appelblad
                       "DE": 28756.6682040598,
                       "PL": 1, "PR": 1, "NL": 2,
                       "mult": 2,
                       "A_values": [1,
                                    0.1,
                                    0.1,
                                    0,
                                    0],}

A2Pi_LSZ = {"RE": 2.1209,
            "BETA": 0.8,
            "GAMMA": 0.1,
            "P": 1,
            "A_values": [1,
                         1],
            "AINF": 1}

A2Pi_bob_rot_params = {
    "RE": 2.1209,
    "RREF": -1.0,
    "P": 4.0,
    "A_values": [0,
     0,
     ],
    "AINF": 0
}


A2Pi_Spin_rot_params = {"RE": 2.1209,
                        "RREF": -1.0,
                        "P": 1.0,
                        "A_values": [0,
                                     0,
                                     ],
                  "AINF": 0,
                  }

A2Pi_lambda_opq_params = {"RE": 2.1209,
                          "P": 1,
                          "A_values": [0,0,0],
                          "AINF": 0}

A2Pi_lambda_p2q_params = {"RE": 2.1209,
                          "P": 1,
                          "A_values": [0,
                                       0,
                                       ],
                          "AINF": 0}

A2Pi_lambda_q_params = {"RE": 2.1209,
                          "P": 1,
                          "A_values": [0,
                                       0,
                                       ],
                          "AINF": 0}


    
    # A2Pi = ElectronicState('A2Pi')
    
    # A2Pi.potential(PotentialEnergy(3,"A2Pi", mult=2, lambda_=1, functional_instance = EMO(**A2Pi_Potential_dict)))
    # A2Pi.potential_energy.set_ab_initio_data(pd.read_csv('/home/shaun/PhD_Work/Linelist_SrH_BaH/Molpro/SrH/SrH_2_2/finalData/SrH_A_Poten.csv'))
    # A2Pi.spinorbitcoupling(SpinOrbitX(functional_instance=POLYNOM_DECAY_24(**A2Pi_LSZ), potential_energy_instance=A2Pi.potential_energy))
    # A2Pi.spin_orbit_coupling.set_ab_initio_data(pd.read_csv('/home/shaun/PhD_Work/Linelist_SrH_BaH/Molpro/SrH/SrH_2_2/finalData/LSZ_A_A.csv'))
    # A2Pi.spinrot(SpinRot(BOBLEROY(**A2Pi_Spin_rot_params),potential_energy_instance=A2Pi.potential_energy))
    # A2Pi.bobrot(BobRot(functional_instance = BOBLEROY(**A2Pi_bob_rot_params), potential_energy_instance=A2Pi.potential_energy))
    
    # A2Pi.lambdap2q(Lambda_p2q(BOBLEROY(**A2Pi_lambda_p2q_params), potential_energy_instance=A2Pi.potential_energy))
    # A2Pi.lambdaq(Lambda_q(BOBLEROY(**A2Pi_lambda_q_params), potential_energy_instance=A2Pi.potential_energy))
    
    # system.add_state(A2Pi)

A2Pi = ElectronicState('A2Pi')

A2Pi.potential(PotentialEnergy(3,"A2Pi", mult=2, lambda_=1, functional_instance = EMO(**prev_A2Pi['poten'])))
A2Pi.potential_energy.set_ab_initio_data(pd.read_csv('./SrH_A_Poten.csv'))
A2Pi.spinorbitcoupling(SpinOrbitX(functional_instance=POLYNOM_DECAY_24(**prev_A2Pi['spin-orbit']), potential_energy_instance=A2Pi.potential_energy))
A2Pi.spin_orbit_coupling.set_ab_initio_data(pd.read_csv('./LSZ_A_A.csv'))
A2Pi.spinrot(SpinRot(POLYNOM(**prev_A2Pi['spin-rot']),potential_energy_instance=A2Pi.potential_energy))
A2Pi.bobrot(BobRot(functional_instance = BOBLEROY(**prev_A2Pi['bob-rot']), potential_energy_instance=A2Pi.potential_energy))

A2Pi.lambdap2q(Lambda_p2q(BOBLEROY(**prev_A2Pi['lambda-p2q']), potential_energy_instance=A2Pi.potential_energy))
A2Pi.lambdaq(Lambda_q(BOBLEROY(**prev_A2Pi['lambda-q']), potential_energy_instance=A2Pi.potential_energy))
A2Pi.lambdaopq(Lambda_opq(BOBLEROY(**A2Pi_lambda_opq_params), potential_energy_instance = A2Pi.potential_energy))

system.add_state(A2Pi)
#%% B2Sigma+ Object


prev_B2Sig = {'poten': prev['poten']['B2Sig+'],
              'bob-rot': prev['bob-rot']['<B2Sig+|BR|B2Sig+>'],
              'spin-rot': prev['spin-rot']['<B2Sig+|SR|B2Sig+>']}

B2Sig_Potential_dict = {"V0": 14312.691,
                       "RE": 2.1169,#86Appelblad
                       "DE": 28756.6682040598,
                       "PL": 1, "PR": 1, "NL": 2,
                       "mult": 2,
                       "A_values": [1,
                                    0.1,
                                    0.1,
                                    0,
                                    0],}

B2Sig_bob_rot_params = {
    "RE": 2.1209,
    "RREF": -1.0,
    "P": 4.0,
    "A_values": [0,
     0,
     ],
    "AINF": 0
}


B2Sig_spin_rot_params = {"RE": 2.1209,
                        "RREF": -1.0,
                        "P": 1.0,
                        "A_values": [0,
                                     0,
                                     ],
                  "AINF": 0,
                  }

B2Sig = ElectronicState('B2Sig+')

# B2Pi.potential(PotentialEnergy(3,"A2Pi", mult=2, lambda_=1, functional_instance = EMO(**prev_B2Pi['poten'])))
# B2Pi.potential_energy.set_ab_initio_data(pd.read_csv('/home/shaun/PhD_Work/Linelist_SrH_BaH/Molpro/SrH/SrH_2_2/finalData/SrH_A_Poten.csv'))

# B2Pi.spinorbitcoupling(SpinOrbitX(functional_instance=POLYNOM_DECAY_24(**prev_B2Pi['spin-orbit']), potential_energy_instance=A2Pi.potential_energy))
# B2Pi.spin_orbit_coupling.set_ab_initio_data(pd.read_csv('/home/shaun/PhD_Work/Linelist_SrH_BaH/Molpro/SrH/SrH_2_2/finalData/LSZ_A_A.csv'))

# B2Pi.spinrot(SpinRot(BOBLEROY(**prev_ABPi['spin-rot']),potential_energy_instance=A2Pi.potential_energy))
# B2Pi.bobrot(BobRot(functional_instance = BOBLEROY(**prev_B2Pi['bob-rot']), potential_energy_instance=A2Pi.potential_energy))

B2Sig.potential(PotentialEnergy(4,"B2Sig+", mult=2, lambda_=0, functional_instance = EMO(**prev_B2Sig['poten']),symmetry="+"))
B2Sig.potential_energy.set_ab_initio_data(pd.read_csv('/home/shaun/pyduoFitting/pyduo-dev/PEC_B2SIG.csv'))

B2Sig.spinrot(SpinRot(POLYNOM(**prev_B2Sig['spin-rot']),potential_energy_instance=B2Sig.potential_energy))
B2Sig.bobrot(BobRot(functional_instance = BOBLEROY(**prev_B2Sig['bob-rot']), potential_energy_instance=B2Sig.potential_energy))

system.add_state(B2Sig)
#%%

C2Sig_Dict_EMO = B2Sig_Potential_dict = {"V0": 21062.4064712496,
                       "RE": 2.72,#86Appelblad
                       "DE": 33157.3081717798,
                       "PL": 10, "PR": 2, "NL": 0,
                       "mult": 2,
                       "A_values": [0.65,
                                    2,
                                    1,
                                    1,
                                    0],}

C2Sig = ElectronicState('C2Sig+')
C2Sig.potential(PotentialEnergy(5,"C2Sig+", mult=2, lambda_=0, functional_instance = EMO(**C2Sig_Dict_EMO),symmetry="+"))
C2Sig.potential_energy.set_ab_initio_data(pd.read_csv('./PEC_C2SIG.csv'),"1e12")
#system.add_state(C2Sig)


#%%

A2Pi_E2Pi_LSZ = {"RE": 2.6,
            "BETA": .8,
            "GAMMA": .02,
            "P": 1,
            "A_values": [1,
                         ],
            "AINF": 1}

X2Pi_E2Pi_LSZ = {"RE": 1.6,
            "BETA": .8,
            "GAMMA": .02,
            "P": 1,
            "A_values": [1,
                         ],
            "AINF": 1}

X2Pi_A2Pi_LSZ = {"RE": 1.8,
            "BETA": .8,
            "GAMMA": 0.02,
            "P": 1,
            "A_values": [1,
                         ],
            "AINF": 1}

B2Pi_A2Pi_LSZ = {"RE": 2.1,
            "BETA": 0.8,
            "GAMMA": 0.02,
            "P": 1,
            "A_values": [1.81914268997923E+00,
                         ],
            "AINF": 1}

B2Pi_E2Pi_LSZ = {"RE": 3.2,
            "BETA": 0.8,
            "GAMMA": 0.02,
            "P": 1,
            "A_values": [-9.59431898900021E-01,
                         ],
            "AINF": 1}

BA_SS = {"RE": 2.102,
            "P": 1,
            "A_values": [0,
                         ],
            "AINF": 0}

BC_Diabatic = {'RE':2.8,
               'GAMMA':0.2,
               'A_values':[0.1]}

off_diag = OffDiagonalComponents()

#off_diag.generate_nac(state1 = B2Sig, state2 = C2Sig, functional_instance = LORENTZ(**BC_Diabatic))

off_diag.generate_ls_x_coupling(state1 = A2Pi, state2 = E2Pi, ab_initio = './LSZ_A_E.csv',
                              functional_instance = POLYNOM_DECAY_24(**A2Pi_E2Pi_LSZ),
                              lambdas = (1,1), sigmas=(0.5,0.5))

off_diag.generate_ls_x_coupling(state1 = electronic_state, state2 = A2Pi, ab_initio = './LSY_A_X.csv',
                              functional_instance = POLYNOM_DECAY_24(**X2Pi_A2Pi_LSZ),
                              lambdas = (0,1), sigmas=(0.5,-0.5))
   
off_diag.generate_ls_x_coupling(state1 = electronic_state, state2 = E2Pi, ab_initio = './LSY_E_X.csv',
                              functional_instance = POLYNOM_DECAY_24(**X2Pi_E2Pi_LSZ),
                              lambdas = (0,1), sigmas=(0.5,-0.5))
   
off_diag.generate_ls_x_coupling(state1 = B2Sig, state2 = E2Pi, ab_initio = './LSY_B_E.csv',
                              functional_instance = POLYNOM_DECAY_24(**B2Pi_E2Pi_LSZ),
                              lambdas = (0,1), sigmas=(0.5,-0.5))

off_diag.generate_ls_x_coupling(state1 = B2Sig, state2 = A2Pi, ab_initio = './LSY_B_A.csv',
                              functional_instance = POLYNOM_DECAY_24(**B2Pi_A2Pi_LSZ),
                              lambdas = (0,1), sigmas=(0.5,-0.5))

# off_diag.generate_spin_spin_couplings(B2Sig, A2Pi, BOBLEROY(**BA_SS), sigmas=(0.5,-0.5))
off_diag.generate_lx(state1=B2Sig, state2=E2Pi, functional_instance=POLYNOM(**{'A_values':[1.11221690742527E+00,3.83739741078587E-02,0]}), lx = (0,'-i'), factor='i')
off_diag.generate_lx(state1=B2Sig, state2=A2Pi, functional_instance=POLYNOM(**{'A_values':[1.15155775139138E+00,0,0]}), lx = (0,'-i'), factor='i')
#off_diag.generate_lx(state1=electronic_state, state2=A2Pi, functional_instance=POLYNOM(**{'A_values':[8.29091849941239E-01,0,0]}), lx = (0,'-i'), factor='i')
#off_diag.generate_lx(state1=electronic_state, state2=E2Pi, functional_instance=POLYNOM(**{'A_values':[0,0,0]}), lx = (0,'-i'), factor='i')
#off_diag.generate_lx(state1=B2Sig, state2=E2Pi, functional_instance=POLYNOM(**{'A_values':[0,0,0]}), lx = (0,'-i'), factor='i')
system.off_diagonals(off_diag)
#%% DUO Setup

system.set_fitting(["0.5 - 27.5"], "SrH_X_1",
                   './XBAE_VIXPGO_41EXP.csv',
                   'frequencies')

system.generate_abinitio_blocks(units="cm-1", fit_factor=1e-12)
#system.run_duo()

#vary_parameters = [f"A{i}" for i in range(6)]
#

xvary= {
        'potential_energy':[f'A{i}' for i in range(0,10)],
        'spin_rot': [f'A{i}' for i in range(4)],
        'bob_rot':[f'A{i}' for i in range(6)]+['AINF'],
       }

bvary= {
        'potential_energy':['V0','DE','RE']+[f'A{i}' for i in range(6)],
        'spin_rot':[f'A{i}' for i in range(7)],
        'bob_rot':[f'A{i}' for i in range(4)]+[f'B{i}' for i in range(4)]+['AINF'],
       }

evary = {#'lambda_opq': ['AINF']+[f'A{i}' for i in range(3)],
         'potential_energy':['V0']+[f'A{i}' for i in range(3)],
         'spin_orbit_coupling':['RE','AINF']+[f'A{i}' for i in range(4)],
         'spin_rot':[f'A{i}' for i in range(5)]+['AINF'],
         'bob_rot': [f'A{i}' for i in range(5)]+['AINF'],
         'lambda_q': ['AINF']+[f'A{i}' for i in range(3)],
         'lambda_p2q': ['AINF']+[f'A{i}' for i in range(3)],
         }

avary = {
        'potential_energy':['V0','RE']+[f'A{i}' for i in range(2)],
        'spin_orbit_coupling':[f'A{i}' for i in range(1)],
        'spin_rot':[f'A{i}' for i in range(4)]+['AINF'],
        'bob_rot':[f'A{i}' for i in range(4)]+['AINF'],
        'lambda_q':[f'A{i}' for i in range(3)]+['AINF'],
        'lambda_p2q':[f'A{i}' for i in range(3)]+['AINF'],    
         }

cvary = {'potential_energy':[f'A{i}' for i in range(5)]}

states_to_vary = {"E2Pi": evary,
                  "X2Sig+": xvary,
                  "A2Pi": avary,
                  "B2Sig+":bvary,

                }

off_diagonal_to_vary = {#'<B2Sig+|DC|C2Sig+>': [f'A{i}' for i in range(2)],
                        
                        #"<X2Sig+|LX|A2Pi>": ['A1','A2'],
                        #"<X2Sig+|LX|E2Pi>": ['A1','A2'],
                        "<B2Sig+|LX|A2Pi>": ['A1','A2'],
                        "<B2Sig+|LX|E2Pi>": ['A1','A2'],
                        "<B2Sig+|LSY|A2Pi>": [f'A{i}' for i in range(4)],
                        '<A2Pi|LSY|E2Pi>': [f'A{i}' for i in range(2)],
                        #'<X2Sig+|LSY|E2Pi>': [f'A{i}' for i in range(2)],
                        #'<X2Sig+|LSY|A2Pi>': [f'A{i}' for i in range(2)],
                        
                        '<B2Sig+|LSY|E2Pi>': [f'A{i}' for i in range(2)],
                        }



#%%

system.set_parameters_to_vary(state_vary = states_to_vary, off_diagonal_vary=off_diagonal_to_vary)
system.run_duo()
# Perform the fitting
system.iterative_fit(state_fit=True,vary_diag_first=True, state_vary = states_to_vary, off_diagonal_vary = off_diagonal_to_vary)

#%%
#pyduo.api.run_duo_input('./SrH_v362.inp')
#%%
pyduo.api.plot_residuals('./SrH_v376.en', energy=True)

#%%
plot_abinitio_calc(pot_file_dataframes('SrH_v376.pot'))

#%%

# en_files = glob.glob('./*.en')
# pngs = glob.glob('./*.png')
# for file in en_files:
#     print(f'Plotting {file}')
#     if file.replace('en','png') not in pngs:
#         try:
#             pyduo.api.plot_residuals(file, energy=True, save=True)
#         except:
#             pass
#         #pyduo.api.run_duo_input(file.replace('en','inp'))
        #pyduo.api.plot_residuals(file, energy=True, save=True)
