#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  6 16:45:06 2023

@author: shaun
"""

class PotentialEnergy:
    def __init__(self, state_number, name, lambda_, mult, functional_instance, sigma=None, symmetry=None):
        self.object_type = "poten"
        self.state_number = state_number
        self.name = name
        self.lambda_ = lambda_
        self.symmetry = symmetry
        self.mult = mult
        self.expansion_type = functional_instance.__class__.__name__
        self.values_dict = functional_instance.to_dict()
        self.ab_initio_data = None  # Placeholder for ab initio data
        self.sigma = sigma

    def set_parameters_to_vary(self, vary_parameters):
        for param in vary_parameters:
            if param in self.values_dict:
                value, _ = self.values_dict[param]
                self.values_dict[param] = (value, True)

    def get_varying_parameters(self):
        """Retrieve parameters that are marked to vary."""
        return {key: val for key, (val, vary) in self.values_dict.items() if vary}

    def update_values(self, updated_values_dict):
        """Update the values after optimization."""
        for key, (value, vary) in self.values_dict.items():
            if vary:
                # Update the value while keeping the vary flag unchanged
                self.values_dict[key] = (updated_values_dict[key], vary)
    
    def output_params_dict(self):
        out_dict = {}
        A_values = []
        for key, value in self.values_dict.items():
            if "A" not in key:
                out_dict[key] = value[0]
            elif "A" in key and key != "AINF":
                A_values.append(value[0])
        out_dict["A_values"] = A_values
        if "AINF" in self.values_dict.keys():
            out_dict['AINF'] = self.values_dict['AINF'][0]
        return out_dict
    
    def set_ab_initio_data(self, df):
        """
        Sets the ab initio data for the potential energy object.
        Expects a pandas DataFrame with columns ['r', 'energy'].
        """
        self.ab_initio_data = df
        
    def create_ab_initio_block(self):
        preamble = f"""abinitio poten {self.state_number}
name {self.name}
lambda {self.lambda_}
type grid
mult {self.mult}
units angstrom cm-1
fit_factor 1e-12
values"""
        values_str = self.ab_initio_data.to_csv(sep='\t', index=False, header=False)
        final_str = '\n'.join([preamble, values_str[:-1],"end"])
        return final_str
        
    def __str__(self):
        lines = []
        lines.append(f"POTEN {self.state_number}")
        lines.append(f'name "{self.name}"')
        lines.append(f"lambda {self.lambda_}")
        
        if self.lambda_ == 0:
            lines.append(f"symmetry {self.symmetry}")
            
        lines.extend([
            f"mult   {self.mult}",
            f"type   {self.expansion_type}",
            "values"
        ])
        
        # Add values from values_dict
        for key, (val, _), in self.values_dict.items():
            lines.append(f"{key}           {val:.14E}")

        lines.append("end")
        return "\n".join(lines)

class SpinOrbitX:
    """
    Spin-orbit coupling object for molpro Brett-Pauli LSZ. Used for intrastate coupling.

    Parameters
    ----------
    functional_instance : TYPE
        DESCRIPTION.
    state_numbers : TYPE, optional
        DESCRIPTION. The default is (0,0).
    name : TYPE, optional
        DESCRIPTION. The default is "<State1|LS|State2>".
    mult : TYPE, optional
        DESCRIPTION. The default is 2.
    lambda_ : TYPE, optional
        DESCRIPTION. The default is (0,0).
    potential_energy_instance : TYPE, optional
        DESCRIPTION. The default is None.
    sigma : TYPE, optional
        DESCRIPTION. The default is (0,0).
    factor : TYPE, optional
        DESCRIPTION. The default is 1.


    """
    def __init__(self, functional_instance, state_numbers=(0,0), name="<State1|LS|State2>",
                 mult=2, lambda_=(0,0), potential_energy_instance=None, sigma=(0,0), factor="i", spin=(0,0)):
        self.object_type = "spin-orbit-x"
        self.potential_energy_instance = potential_energy_instance
        self.functional_instance = functional_instance
        self.expansion_type = functional_instance.__class__.__name__
        self.state_numbers = state_numbers
        self.sigma = sigma
        self.factor = factor
        if potential_energy_instance:
            self.name = f'<{self.potential_energy_instance.name}|LS|{self.potential_energy_instance.name}>'
            self.lambda_ = (self.potential_energy_instance.lambda_, self.potential_energy_instance.lambda_)
            self.mult = self.potential_energy_instance.mult
            self.state_numbers = (self.potential_energy_instance.state_number,self.potential_energy_instance.state_number)
            self.spin = ((self.mult-1)/2,(self.mult-1)/2)
            self.sigma = self.spin
        else:
            self.name = name
            self.mult = mult
            self.lambda_ = lambda_
            self.spin = spin
            self.sigma = sigma
        self.values_dict = self.functional_instance.to_dict()
    
    def set_ab_initio_data(self, df):
        self.ab_initio_data = df
    
    def create_ab_initio_block(self):
        preamble = f"""abinitio spin-orbit-x {self.state_numbers[0]} {self.state_numbers[1]}
name {self.name}
spin {self.spin[0]} {self.spin[0]}
lambda {self.lambda_[0]} {self.lambda_[1]}
sigma {self.sigma} {self.sigma}
factor i
<x|LZ|y> -i -i
type grid
units angstrom cm-1
fit_factor 1-e12
values"""
        values_str = self.ab_initio_data.to_csv(sep='\t', index=False, header=False)
        final_str = '\n'.join([preamble, values_str[:-1],"end"])
        return final_str
    
    def set_parameters_to_vary(self, vary_parameters):
        for param in vary_parameters:
            if param in self.values_dict:
                value, _ = self.values_dict[param]
                self.values_dict[param] = (value, True)
    
    def get_varying_parameters(self):
        """Retrieve parameters that are marked to vary."""
        return {key: val for key, (val, vary) in self.values_dict.items() if vary}

    def update_values(self, updated_values_dict):
        """Update the values after optimization."""
        for key, (value, vary) in self.values_dict.items():
            if vary:
                # Update the value while keeping the vary flag unchanged
                self.values_dict[key] = (updated_values_dict[key], vary)
    
    def output_params_dict(self):
        out_dict = {}
        A_values = []
        for key, value in self.values_dict.items():
            if "A" not in key:
                out_dict[key] = value[0]
            elif "A" in key and key != "AINF":
                A_values.append(value[0])
        out_dict["A_values"] = A_values
        if "AINF" in self.values_dict.keys():
            out_dict['AINF'] = self.values_dict['AINF'][0]
        return out_dict            

    def __str__(self):
        lines = []
        #values_str = "\n".join([f'{key:<13}{value[0]:.14E}' for key, value in self.values_dict.items()])
        preamble =  f"""spin-orbit-x {self.state_numbers[0]} {self.state_numbers[1]}
name "{self.name}"
spin {self.spin[0]} {self.spin[0]}
lambda {self.lambda_[0]} {self.lambda_[0]}
sigma {self.sigma[0]} {self.sigma[1]}
factor {self.factor}
<x|LZ|y> -i -i
type {self.expansion_type}
units angstrom cm-1"""
        lines.append(preamble)
        if self.expansion_type == "POLYNOM_DECAY_24":
            lines.append('morphing')
        lines.append('values')
        for key, (val, _) in self.values_dict.items():
            lines.append(f"{key}           {val:.14E}")
        lines.append('end')
        return "\n".join(lines)
    

class SpinOrbit(SpinOrbitX):
    def __init__(self, functional_instance, spin=None, sigma=None, lambda_=None,
                 name="<State1|LS|State2>",
                 potential_energy_instance1=None, potential_energy_instance2=None,
                 factor="i"):
        self.object_type = "spin-orbit"
        self.potential_energy_instance1 = potential_energy_instance1
        self.potential_energy_instance2 = potential_energy_instance2
        self.functional_instance = functional_instance
        self.expansion_type = functional_instance.__class__.__name__
        self.factor = factor
        if potential_energy_instance1 and potential_energy_instance2:
            self.name = f'<{self.potential_energy_instance1.name}|LSY|{self.potential_energy_instance2.name}>'
            self.mult = (self.potential_energy_instance1.mult, self.potential_energy_instance2.mult)
            self.state_numbers = (self.potential_energy_instance1.state_number,self.potential_energy_instance2.state_number)
        else:
            raise ValueError("Requires both potential_energy_instance1 and potential_energy_instance2")
        self.values_dict = self.functional_instance.to_dict()
        if not spin:
            self.spin = ((self.mult[0]-1)/2,(self.mult[1]-1)/2)
        else:
            self.spin = spin
        if not lambda_:
            self.lambda_ = (self.potential_energy_instance1.lambda_, self.potential_energy_instance2.lambda_)
        else:
            self.lambda_ = lambda_
        if not sigma:
            self.sigma = (self.spin[0], -self.spin[1])
        else:
            self.sigma = sigma
        
    def create_ab_initio_block(self):
        preamble = f"""abinitio spin-orbit {self.state_numbers[0]} {self.state_numbers[1]}
name {self.name}
spin {self.spin[0]} {self.spin[0]}
lambda {self.lambda_[0]} {self.lambda_[1]}
sigma {self.sigma} {self.sigma}
factor i
type grid
units angstrom cm-1
fit_factor 1-e12
values"""
        values_str = self.ab_initio_data.to_csv(sep='\t', index=False, header=False)
        final_str = '\n'.join([preamble, values_str[:-1],"end"])
        return final_str    
    
    def __str__(self):
        lines = []
        #values_str = "\n".join([f'{key:<13}{value[0]:.14E}' for key, value in self.values_dict.items()])
        preamble =  f"""spin-orbit {self.state_numbers[0]} {self.state_numbers[1]}
name "{self.name}"
spin {self.spin[0]} {self.spin[1]}
lambda {self.lambda_[0]} {self.lambda_[1]}
sigma {self.sigma[0]} {self.sigma[1]}
factor {self.factor}
type {self.expansion_type}
units angstrom cm-1"""
        lines.append(preamble)
        if self.expansion_type == "POLYNOM_DECAY_24":
            lines.append('morphing')
        lines.append('values')
        for key, (val, _) in self.values_dict.items():
            lines.append(f"{key}           {val:.14E}")
        lines.append('end')
        return "\n".join(lines)

    
class BobRot:
    def __init__(self, functional_instance, state_numbers=(0, 0), name="<State1|BR|State2>",
                 mult=2, lambda_=(0, 0),potential_energy_instance=None):
        self.potential_energy_instance = potential_energy_instance
        self.functional_instance = functional_instance
        self.expansion_type = functional_instance.__class__.__name__
        self.state_numbers = state_numbers
        if potential_energy_instance:
            self.name = f'<{self.potential_energy_instance.name}|BR|{self.potential_energy_instance.name}>'
            self.lambda_ = (self.potential_energy_instance.lambda_, self.potential_energy_instance.lambda_)
            self.mult = self.potential_energy_instance.mult
            self.state_numbers = (self.potential_energy_instance.state_number,self.potential_energy_instance.state_number)
        else:
            self.name = name
            self.mult = mult
            self.lambda_ = lambda_
        self.values_dict = self.functional_instance.to_dict()
    
    def set_parameters_to_vary(self, vary_parameters):
        for param in vary_parameters:
            if param in self.values_dict:
                value, _ = self.values_dict[param]
                self.values_dict[param] = (value, True)
    
    def get_varying_parameters(self):
        """Retrieve parameters that are marked to vary."""
        return {key: val for key, (val, vary) in self.values_dict.items() if vary}

    def update_values(self, updated_values_dict):
        """Update the values after optimization."""
        for key, (value, vary) in self.values_dict.items():
            if vary:
                # Update the value while keeping the vary flag unchanged
                self.values_dict[key] = (updated_values_dict[key], vary)
     
    def output_params_dict(self):
        out_dict = {}
        A_values = []
        for key, value in self.values_dict.items():
            if "A" not in key:
                out_dict[key] = value[0]
            elif "A" in key and key != "AINF":
                A_values.append(value[0])
        out_dict["A_values"] = A_values
        if "AINF" in self.values_dict.keys():
            out_dict['AINF'] = self.values_dict['AINF'][0]
        return out_dict           
     
    def __str__(self):
        values_str = "\n".join([f"{key:<13}{value[0]:.14E}" for key, value in self.values_dict.items()])
        lines = []
        lines.append(f'Bob-Rot {self.state_numbers[0]} {self.state_numbers[1]}')
        lines.append(f'name {self.name}')
        lines.append(f'mult {self.mult}')
        lines.append(f'lambda {self.lambda_[0]} {self.lambda_[1]}')
        lines.append(f'type {self.expansion_type}')
        lines.append(f'factor 1.0')
        lines.append(f'values')
        for key, (val, _) in self.values_dict.items():
            lines.append(f"{key}           {val:.14E}")
        lines.append("end")
        return "\n".join(lines)
    
class SpinRot:
    def __init__(self, functional_instance, state_numbers=(0,0), name="<State1|SR|State2>",
                 lambda_=(0,0), spin=(0,0), factor = 1.0, units = 'cm-1', potential_energy_instance=None):
        
        self.potential_energy_instance = potential_energy_instance
        self.functional_instance = functional_instance
        self.expansion_type = functional_instance.__class__.__name__
        self.state_numbers = state_numbers
        self.factor=factor
        self.units = 'cm-1'
        if potential_energy_instance:
            self.name = f'<{self.potential_energy_instance.name}|SR|{self.potential_energy_instance.name}>'
            self.lambda_ = (self.potential_energy_instance.lambda_, self.potential_energy_instance.lambda_)
            self.state_numbers = (self.potential_energy_instance.state_number,self.potential_energy_instance.state_number)
            self.mult = self.potential_energy_instance.mult
            self.spin = ((self.mult-1)/2,(self.mult-1)/2)
        else:
            self.spin=spin
            self.name = name
            self.lambda_ = lambda_
        self.values_dict = self.functional_instance.to_dict()
        
        
    def set_parameters_to_vary(self, vary_parameters):
        for param in vary_parameters:
            if param in self.values_dict:
                value, _ = self.values_dict[param]
                self.values_dict[param] = (value, True)
    
    def get_varying_parameters(self):
        """Retrieve parameters that are marked to vary."""
        return {key: val for key, (val, vary) in self.values_dict.items() if vary}

    def update_values(self, updated_values_dict):
        """Update the values after optimization."""
        for key, (value, vary) in self.values_dict.items():
            if vary:
                # Update the value while keeping the vary flag unchanged
                self.values_dict[key] = (updated_values_dict[key], vary)
     
    def output_params_dict(self):
        out_dict = {}
        A_values = []
        for key, value in self.values_dict.items():
            if "A" not in key:
                out_dict[key] = value[0]
            elif "A" in key and key != "AINF":
                A_values.append(value[0])
        out_dict["A_values"] = A_values
        if "AINF" in self.values_dict.keys():
            out_dict['AINF'] = self.values_dict['AINF'][0]
        return out_dict           
    
    def __str__(self):
        values_str = "\n".join([f"{key:<13}{value[0]:.14E}" for key, value in self.values_dict.items()])
        lines = []
        lines.append(f'spin-rot {self.state_numbers[0]} {self.state_numbers[1]}')
        lines.append(f'name {self.name}')
        lines.append(f'spin {self.spin[0]} {self.spin[1]}')
        lines.append(f'lambda {self.lambda_[0]} {self.lambda_[1]}')
        lines.append(f'factor {self.factor}')
        lines.append(f'type {self.expansion_type}')
        lines.append(f'values')
        for key, (val, _) in self.values_dict.items():
            lines.append(f"{key}           {val:.14E}")
        lines.append("end")
        return "\n".join(lines)
    
class SpinSpin(SpinRot):
    def __init__(self):
        self.type="spin_rot"
    def set_parameters_to_vary(self, vary_parameters):
        for param in vary_parameters:
            if param in self.values_dict:
                value, _ = self.values_dict[param]
                self.values_dict[param] = (value, True)
    
    def get_varying_parameters(self):
        """Retrieve parameters that are marked to vary."""
        return {key: val for key, (val, vary) in self.values_dict.items() if vary}

    def update_values(self, updated_values_dict):
        """Update the values after optimization."""
        for key, (value, vary) in self.values_dict.items():
            if vary:
                # Update the value while keeping the vary flag unchanged
                self.values_dict[key] = (updated_values_dict[key], vary)
                
    def output_params_dict(self):
        out_dict = {}
        A_values = []
        for key, value in self.values_dict.items():
            if "A" not in key:
                out_dict[key] = value[0]
            elif "A" in key and key != "AINF":
                A_values.append(value[0])
        out_dict["A_values"] = A_values
        if "AINF" in self.values_dict.keys():
            out_dict['AINF'] = self.values_dict['AINF'][0]
        return out_dict        
        
    def __str__(self):
        values_str = "\n".join([f"{key:<13}{value[0]:.14E}" for key, value in self.values_dict.items()])
        lines = []
        lines.append(f'spin-spin {self.state_numbers[0]} {self.state_numbers[1]}')
        lines.append(f'name {self.name}')
        lines.append(f'spin {self.spin[0]} {self.spin[1]}')
        lines.append(f'lambda {self.lambda_[0]} {self.lambda_[1]}')
        lines.append(f'sigma {(self.mult-1)/2} {(self.mult-1)/2}')
        lines.append(f'factor {self.factor}')
        lines.append(f'type {self.expansion_type}')
        lines.append(f'values')
        for key, (val, _) in self.values_dict.items():
            lines.append(f"{key}           {val:.14E}")
        lines.append("end")
        return "\n".join(lines)
    
class Lambda_opq:
    def __init__(self, functional_instance, state_numbers = (0,0), name = "<State1|lambda-opq|State2>",
                 lambda_ = (0,0), spin = (0,0), sigma = (-0.5, 0.5), factor = 1, potential_energy_instance = None):
        self.potential_energy_instance = potential_energy_instance
        self.functional_instance = functional_instance
        self.expansion_type = functional_instance.__class__.__name__
        self.state_numbers = state_numbers
        self.factor = factor
        self.sigma = sigma
        if potential_energy_instance:
            self.name = f'<{self.potential_energy_instance.name}|lambda-opq|{self.potential_energy_instance.name}>'
            self.lambda_ = (self.potential_energy_instance.lambda_, self.potential_energy_instance.lambda_)
            self.state_numbers = (self.potential_energy_instance.state_number,self.potential_energy_instance.state_number)
            self.spin = (self.potential_energy_instance.mult-1)/2
            self.sigma = (-self.spin, self.spin)
        else:
            self.name = name
            self.lambda_ = lambda_
            self.sigma = sigma
        self.values_dict = self.functional_instance.to_dict()
        
    def set_parameters_to_vary(self, vary_parameters):
        for param in vary_parameters:
            if param in self.values_dict:
                value, _ = self.values_dict[param]
                self.values_dict[param] = (value, True)
    
    def get_varying_parameters(self):
        """Retrieve parameters that are marked to vary."""
        return {key: val for key, (val, vary) in self.values_dict.items() if vary}

    def update_values(self, updated_values_dict):
        """Update the values after optimization."""
        for key, (value, vary) in self.values_dict.items():
            if vary:
                # Update the value while keeping the vary flag unchanged
                self.values_dict[key] = (updated_values_dict[key], vary)    
    
    def output_params_dict(self):
        out_dict = {}
        A_values = []
        for key, value in self.values_dict.items():
            if "A" not in key:
                out_dict[key] = value[0]
            elif "A" in key and key != "AINF":
                A_values.append(value[0])
        out_dict["A_values"] = A_values
        if "AINF" in self.values_dict.keys():
            out_dict['AINF'] = self.values_dict['AINF'][0]
        return out_dict
    
    def __str__(self):
        lines = []
        values_str = "\n".join([f'{key:<13}{value[0]:.14E}' for key, value in self.values_dict.items()])
        preamble = f"""lambda-opq {self.state_numbers[0]} {self.state_numbers[1]}
name "{self.name}"
spin {self.spin} {self.spin}
sigma {self.sigma[0]} {self.sigma[1]}
type {self.expansion_type}
units cm-1
lambda {self.lambda_[0]} {self.lambda_[1]}
values"""
        lines.append(preamble)
        for key, (val, _) in self.values_dict.items():
            lines.append(f"{key}           {val:.14E}")
        lines.append("end")
        return '\n'.join(lines)

class Lambda_p2q:
    def __init__(self, functional_instance, state_numbers = (0,0), name = "<State1|lambda-p2q|State2>",
                 lambda_ = (0,0), spin = (0,0), sigma = (0,0), factor = 1, potential_energy_instance = None):
        self.potential_energy_instance = potential_energy_instance
        self.functional_instance = functional_instance
        self.expansion_type = functional_instance.__class__.__name__
        self.state_numbers = state_numbers
        self.factor = factor
        if potential_energy_instance:
            self.name = f'<{self.potential_energy_instance.name}|lambda-p2q|{self.potential_energy_instance.name}>'
            self.lambda_ = (self.potential_energy_instance.lambda_, self.potential_energy_instance.lambda_)
            self.state_numbers = (self.potential_energy_instance.state_number,self.potential_energy_instance.state_number)
            self.spin = (self.potential_energy_instance.mult-1)/2
            self.sigma = (-self.spin, self.spin)
        else:
            self.name = name
            self.lambda_ = lambda_
            self.sigma = sigma
        self.values_dict = self.functional_instance.to_dict()
    
    def set_parameters_to_vary(self, vary_parameters):
        for param in vary_parameters:
            if param in self.values_dict:
                value, _ = self.values_dict[param]
                self.values_dict[param] = (value, True)
    
    def get_varying_parameters(self):
        """Retrieve parameters that are marked to vary."""
        return {key: val for key, (val, vary) in self.values_dict.items() if vary}

    def update_values(self, updated_values_dict):
        """Update the values after optimization."""
        for key, (value, vary) in self.values_dict.items():
            if vary:
                # Update the value while keeping the vary flag unchanged
                self.values_dict[key] = (updated_values_dict[key], vary)
    
    def output_params_dict(self):
        out_dict = {}
        A_values = []
        for key, value in self.values_dict.items():
            if "A" not in key:
                out_dict[key] = value[0]
            elif "A" in key and key != "AINF":
                A_values.append(value[0])
        out_dict["A_values"] = A_values
        if "AINF" in self.values_dict.keys():
            out_dict['AINF'] = self.values_dict['AINF'][0]
        return out_dict
    
    def __str__(self):
        lines = []
        values_str = "\n".join([f'{key:<13}{value[0]:.14E}' for key, value in self.values_dict.items()])
        preamble = f"""lambda-p2q {self.state_numbers[0]} {self.state_numbers[1]}
name "{self.name}"
spin {self.spin} {self.spin}
sigma {self.sigma[0]} {self.sigma[1]}
type {self.expansion_type}
units cm-1
lambda {self.lambda_[0]} {self.lambda_[1]}
values"""
        lines.append(preamble)
        for key, (val, _) in self.values_dict.items():
            lines.append(f"{key}           {val:.14E}")
        lines.append("end")
        return '\n'.join(lines)

class Lambda_q:
    def __init__(self, functional_instance, state_numbers = (0,0), name = "<State1|lambda-q|State2>",
                 lambda_ = (0,0), spin = (0,0), factor = 1, potential_energy_instance = None):
        self.potential_energy_instance = potential_energy_instance
        self.functional_instance = functional_instance
        self.expansion_type = functional_instance.__class__.__name__
        self.state_numbers = state_numbers
        self.factor = factor
        if potential_energy_instance:
            self.name = f'<{self.potential_energy_instance.name}|lambda-q|{self.potential_energy_instance.name}>'
            self.lambda_ = (self.potential_energy_instance.lambda_, self.potential_energy_instance.lambda_)
            self.state_numbers = (self.potential_energy_instance.state_number,self.potential_energy_instance.state_number)
            self.spin = (self.potential_energy_instance.mult-1)/2
        else:
            self.name = name
            self.lambda_ = lambda_
        self.values_dict = self.functional_instance.to_dict()
        
    def set_parameters_to_vary(self, vary_parameters):
        for param in vary_parameters:
            if param in self.values_dict:
                value, _ = self.values_dict[param]
                self.values_dict[param] = (value, True)
    
    def get_varying_parameters(self):
        """Retrieve parameters that are marked to vary."""
        return {key: val for key, (val, vary) in self.values_dict.items() if vary}

    def update_values(self, updated_values_dict):
        """Update the values after optimization."""
        for key, (value, vary) in self.values_dict.items():
            if vary:
                # Update the value while keeping the vary flag unchanged
                self.values_dict[key] = (updated_values_dict[key], vary)
    
    def output_params_dict(self):
        out_dict = {}
        A_values = []
        for key, value in self.values_dict.items():
            if "A" not in key:
                out_dict[key] = value[0]
            elif "A" in key and key != "AINF":
                A_values.append(value[0])
        out_dict["A_values"] = A_values
        if "AINF" in self.values_dict.keys():
            out_dict['AINF'] = self.values_dict['AINF'][0]
        return out_dict
    
    def __str__(self):
        lines = []
        values_str = "\n".join([f'{key:<13}{value[0]:.14E}' for key, value in self.values_dict.items()])
        preamble = f"""lambda-q {self.state_numbers[0]} {self.state_numbers[1]}
name "{self.name}"
spin {self.spin} {self.spin}
type {self.expansion_type}
factor {self.factor}
units cm-1
lambda {self.lambda_[0]} {self.lambda_[1]}
values"""
        lines.append(preamble)
        for key, (val, _) in self.values_dict.items():
            lines.append(f"{key}           {val:.14E}")
        lines.append("end")
        return '\n'.join(lines)



