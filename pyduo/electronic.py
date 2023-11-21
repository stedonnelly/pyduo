from .containers import *
from .functionals import *
from .components import *
import pandas as pd

class ElectronicState:
    def __init__(self, name):
        self.name = name
        self.potential_energy = None
        self.spin_orbit_coupling = None
        self.spin_rot = None
        self.bob_rot = None
        self.dipole_moments = None
        self.spin_spin = None
        self.lambda_opq = None
        self.lambda_p2q = None
        self.lambda_q = None
        self.off_diagonal = None
        

    def potential(self, potential_energy_object):
        """Sets the potential energy component of the electronic state."""
        self.potential_energy = potential_energy_object
        
    def bobrot(self, bobrot_object):
        self.bob_rot = bobrot_object
    
    def spinrot(self, spinrot_object):
        self.spin_rot = spinrot_object
    # We can add similar methods for other components as we define them
    
    def spinspin(self, spinspin_object):
        self.spin_spin = spinspin_object
        
    def spinorbitcoupling(self, spinorbitcoupling_object):
        self.spin_orbit_coupling = spinorbitcoupling_object
    
    def lambdaopq(self, lambda_opq_object):
        self.lambda_opq = lambda_opq_object
        
    def lambdap2q(self, lambda_p2q_object):
        self.lambda_p2q = lambda_p2q_object
        
    def lambdaq(self, lambda_q_object):
        self.lambda_q = lambda_q_object
    
    def __str__(self):
        self.components = []
        if self.potential_energy:
            self.components.append(str(self.potential_energy))
        if self.spin_orbit_coupling:
            self.components.append(str(self.spin_orbit_coupling))
        if self.bob_rot:
            self.components.append(str(self.bob_rot))
        if self.spin_rot:
            self.components.append(str(self.spin_rot))
        if self.spin_spin:
            self.components.append(str(self.spin_spin))
        if self.lambda_opq:
            self.components.append(str(self.lambda_opq))
        if self.lambda_p2q:
            self.components.append(str(self.lambda_p2q))
        if self.lambda_q:
            self.components.append(str(self.lambda_q))
        if self.off_diagonal:
            for component in self.off_diagonal.components:
                self.components.append(str(component))
        #We'll add the string representations of other components here as they're defined
        return "\n\n".join(self.components)

class SigmaState(ElectronicState):
    def __init__(self, name: str, state_number: int, mult:int, functional_instance=None):
        self.name = name
        self.potential_energy = None
        self.spin_orbit_coupling = None
        self.spin_rot = None
        self.bob_rot = None
        self.dipole_moments = None
        self.spin_spin = None
        self.lambda_opq = None
        self.lambda_p2q = None
        self.lambda_q = None
        self.off_diagonal = None
        if not functional_instance:
            self.functional_details = {"V0": 0,
                                   "RE": 1.4,
                                   "DE": 10000,
                                   "PL": 1, "PR": 1, "NL": 2,
                                   "mult": 2,
                                   "A_values": [1,
                                                0,
                                                0,
                                                0,
                                                0],}
            self.functional_instance = EMO()
        else:
            self.functional_instance = functional_instance
        
        self.potential_energy = PotentialEnergy(state_number, name, 0, mult, self.functional_instance)

class PiState(ElectronicState):
    def __init__(self, name: str, state_number: int, mult:int, functional_instance=None):
        self.name = name
        self.potential_energy = None
        self.spin_orbit_coupling = None
        self.spin_rot = None
        self.bob_rot = None
        self.dipole_moments = None
        self.spin_spin = None
        self.lambda_opq = None
        self.lambda_p2q = None
        self.lambda_q = None
        self.off_diagonal = None
        if not functional_instance:
            self.functional_details = {"V0": 0,
                                   "RE": 1.4,
                                   "DE": 10000,
                                   "PL": 1, "PR": 1, "NL": 2,
                                   "mult": 2,
                                   "A_values": [1,
                                                0,
                                                0,
                                                0,
                                                0],}
            self.functional_instance = EMO()
        else:
            self.functional_instance = functional_instance
        
        self.potential_energy = PotentialEnergy(state_number, name, 1, mult, self.functional_instance)



class OffDiagonalComponents:
    def __init__(self):
        self.components = {}
    
    def generate_ls_coupling(self, state1, state2, ab_initio, functional_instance,
                             sigmas = (0,0),lambdas = (0,0), spin = None,
                             override_ab_initio = False):
        spin_orbit_object = SpinOrbit(functional_instance,potential_energy_instance1 = state1.potential_energy, potential_energy_instance2 = state2.potential_energy,
                                      sigma = sigmas, lambda_ = lambdas, spin = spin)
        if ab_initio:
            spin_orbit_object.set_ab_initio_data(pd.read_csv(ab_initio))
        else:
            if override_ab_initio:
                pass
            else:
                raise ValueError("No ab initio data provided")
        self.components[spin_orbit_object.name] = spin_orbit_object
    
    def generate_ls_x_coupling(self, state1, state2, ab_initio, functional_instance,
                             sigmas = (0,0),lambdas = (0,0), spin = None,
                             override_ab_initio = False):
        spin_orbit_object = SpinOrbitX_Y(functional_instance,potential_energy_instance1 = state1.potential_energy, potential_energy_instance2 = state2.potential_energy,
                                      sigma = sigmas, lambda_ = lambdas, spin = spin)
        if ab_initio:
            spin_orbit_object.set_ab_initio_data(pd.read_csv(ab_initio))
        else:
            if override_ab_initio:
                pass
            else:
                raise ValueError("No ab initio data provided")
        self.components[spin_orbit_object.name] = spin_orbit_object
    
    def generate_spin_spin_couplings(self, state1, state2, functional_instance,sigmas=None):
        spin_spin_coupling = SpinSpin(state1, state2, functional_instance,sigmas)
        self.components[spin_spin_coupling.name] = spin_spin_coupling
        
    def generate_diabatic(self, state1, state2, functional_instance):
        diabatic_object = diabatic(state1, state2,functional_instance)
        self.components[diabatic_object.name] = diabatic_object
    
    def generate_nac(self, state1, state2, functional_instance):
        nac_object = NAC(state1, state2,functional_instance)
        self.components[nac_object.name] = nac_object
        



