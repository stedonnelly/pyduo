#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  2 20:27:20 2023

@author: shaun
"""
from functionals import *
from containers import *
from electronic import *
import os
import subprocess
from api import *
from scipy.optimize import minimize

class MolecularSystem:
    def __init__(self):
        self.states = {}
        self.abinitio = []
        self.molecule_info = None
        self.grid = None
        self.solution_method = None
        self.abinitio_blocks = []
        self.fitting_active = False
        self.diagonalizer = None
        self.off_diagonal = {}
        self.components = []
        
    def add_state(self, electronic_state):
        """Adds an ElectronicState object to the states list."""
        self.states[electronic_state.name] = electronic_state
        
    def off_diagonals(self, off_diagonal_object):
        for key, component in off_diagonal_object.components.items():
            self.off_diagonal[key] = component
        
    def remove_state(self, state_number):
        """Removes an ElectronicState object based on its state number."""
        self.states = [state for state in self.states if state.potential_energy.state_number != state_number]
        
    def get_state(self, state_number):
        """Retrieves an ElectronicState object based on its state number."""
        for state in self.states:
            if state.potential_energy.state_number == state_number:
                return state
        return None
    
    def molecule(self, atom1, atom2, mass1, mass2):
        """Sets the molecule information."""
        self.molecule_info = MoleculeInformation(atom1, atom2, mass1, mass2)
        
    def integration_grid(self, grid_block):
        """Sets the integration grid."""
        self.grid = grid_block
        
    def set_solution_method(self, method):
        """Sets the solution method."""
        self.solution_method = method
    
    def set_fitting(self, jlist, output_name, csv_path,input_type='energies'):
        """Sets the fitting block."""
        input_name, output_name = self.generate_unique_filename()
        if input_type=='energies':
            self.fitting = Fitting(jlist, output_name.split('.')[-1], csv_path)
        if input_type=='frequencies':
            self.fitting = FreqFitting(jlist, output_name.split('.')[-1], csv_path)
        self.weights = self.fitting.df

    def generate_abinitio_blocks(self, units, fit_factor):
        """
        Generates ab initio blocks for each state object and stores them.
        Assumes that ab initio data has been set on each potential energy object.
        """
        for name, state in self.states.items():
            if state.potential_energy and state.potential_energy.ab_initio_data is not None:
                abinitio_block = AbInitioBlock(state.potential_energy, units, fit_factor)
                # Assuming AbInitioBlock can handle a DataFrame directly now:
                abinitio_block.set_values_from_dataframe(state.potential_energy.ab_initio_data)
                self.abinitio_blocks.append(abinitio_block.generate_abinitio())
                
            if state.spin_orbit_coupling and state.spin_orbit_coupling.ab_initio_data is not None:
                abinitio_block = AbInitioBlock(state.spin_orbit_coupling, units, fit_factor)
                # Assuming AbInitioBlock can handle a DataFrame directly now:
                abinitio_block.set_values_from_dataframe(state.spin_orbit_coupling.ab_initio_data)
                self.abinitio_blocks.append(abinitio_block.generate_abinitio())
                
        if self.off_diagonal:
            for key, component in self.off_diagonal.items():
                self.abinitio_blocks.append(component.create_ab_initio_block())

    def generate_unique_filename(self, base_path='.'):
        """
        Generates a unique filename for the input and output files based on the molecule name.
        """
        # Assume the molecule name is stored in self.molecule_info and is of format "Atom1Atom2"
        molecule_name = f"{self.molecule_info.atom1}{self.molecule_info.atom2}"
        version_number = 1
        while True:
            input_filename = f"{base_path}/{molecule_name}_v{version_number}.inp"
            output_filename = f"{base_path}/out.{molecule_name}_v{version_number}"
            # Check if either the input or output file already exists
            if not os.path.exists(input_filename) and not os.path.exists(output_filename):
                return input_filename, output_filename
            else:
                version_number += 1
    def run_duo(self):
        """
        Runs DUO with the generated input file and saves the output.
        """
        input_str = self.generate_input()

    
        # Write the input string to the input file
        with open(self.input_filename, 'w') as f:
            f.write(input_str)
    
        # Execute DUO using the system terminal
        command = f"duo < {self.input_filename} > {self.output_filename}"
        subprocess.run(command, shell=True)
        if "frequencies" in str(self.fitting):
            self.output = duo_output_object(self.output_filename.split('.')[-1]+'.en', weights=self.weights,output_type='frequencies')
        elif "energies" in str(self.fitting):
            self.output = duo_output_object(self.output_filename.split('.')[-1]+'.en', weights=self.weights,output_type='energies')
    
    def generate_input(self):
        """Constructs the final input for DUO by iterating over all stored components."""
        self.components = []
        self.states = dict(sorted(list(self.states.items()), key=lambda x: x[1].potential_energy.state_number))
        # Molecule Information
        if self.molecule_info:
            self.components.append(str(self.molecule_info))
            
        # Total Number of States
        self.components.append(f"nstates {len(self.states)}")
        self.components.append('jrot 0.5 - 1.5\n\nsymmetry Cs(M)')
        
        # Solution Methods
        if self.solution_method:
            self.components.append(f"SolutionMethod {self.solution_method}")
        
        # Integration Grid
        if self.grid:
            self.components.append(str(self.grid))
        
        if self.diagonalizer:
            self.components.append(str(self.diagonalizer))
        # Contraction Block
        potential_blocks = [self.states[state].potential_energy for state in self.states]
        contraction_block = ContractionBlock(potential_blocks)
        self.components.append(str(contraction_block))
        
        # Electronic States
        states_str = [str(state) for name, state in self.states.items()]
        self.components.extend(states_str)
        
        if self.off_diagonal:
            off_diagonals_str = [str(component) for key, component in self.off_diagonal.items()]
            self.components.extend(off_diagonals_str)
        
        if self.fitting:
            self.components.append(str(self.fitting))
        
        # Ab Initio Blocks
        self.components.extend(self.abinitio_blocks)
        
        # We'll expand for Fitting Block and Ab Initio Blocks as we define those components
        
        return "\n\n".join(self.components)


    def extract_parameters(self,state ,component_name):
        """
        Extracts parameters marked as True (for fitting) from the values_dict of the specified component.
        
        :param component_name: Name of the component attribute from which to extract parameters.
        """
        params = []
        # Use getattr to dynamically access the component from the state
        component = getattr(state, component_name, None)
        if component is not None:
            for param, (value, fit) in component.values_dict.items():
                if fit:
                    params.append(value)
        return params
    
    def extract_off_diag_parameters(self,component):
        """
        Extracts parameters marked as True (for fitting) from the values_dict of the specified component.
        
        :param component_name: Name of the component attribute from which to extract parameters.
        """
        params = []
        # Use getattr to dynamically access the component from the state
        if component is not None:
            for param, (value, fit) in component.values_dict.items():
                if fit:
                    params.append(value)
        return params
    
    def update_parameters(self, state, params, component_name):
        """
        Updates the values_dict with a new set of parameters for the specified component.
        Assumes params is a flat list of values corresponding to the order in which
        extract_parameters() extracts the parameters.
        """
        param_index = 0
        component = getattr(state, component_name)  # Access the component by name
        for param, (value, fit) in component.values_dict.items():
            if fit:
                component.values_dict[param] = (params[param_index], fit)
                param_index += 1
    
    def update_off_diag(self, params, component):
        """
        Updates the values_dict with a new set of parameters for the specified component.
        Assumes params is a flat list of values corresponding to the order in which
        extract_parameters() extracts the parameters.
        """
        param_index = 0  # Access the component by name
        for param, (value, fit) in component.values_dict.items():
            if fit:
                component.values_dict[param] = (params[param_index], fit)
                param_index += 1
    
    def set_parameters_to_vary(self, state_vary=None, off_diagonal_vary=None):
        #state.potential_energy.set_parameters_to_vary(vary_parameters) 
        if state_vary:
            for name, state in self.states.items():
                if name in state_vary:
                    for component_name, params in state_vary[name].items():
                        component = getattr(state, component_name, None)
                        
                        # If the component exists and has a method to set parameters to vary
                        if component and hasattr(component, 'set_parameters_to_vary'):
                            # Call the method on the component with the list of parameters to vary
                            component.set_parameters_to_vary(params)
        if off_diagonal_vary:
            for name, off_diag in self.off_diagonal.items():
                if name in off_diagonal_vary:
                    for component_name, params in off_diagonal_vary.items():
                        component = self.off_diagonal[component_name]
                        if component and hasattr(component, 'set_parameters_to_vary'):
                            # Call the method on the component with the list of parameters to vary
                            component.set_parameters_to_vary(params)
    
    def initialise(self):
        self.input_filename, self.output_filename = self.generate_unique_filename()
    
    def single_fit_step(self, state, component_name, verbose=False):
        """
        Performs a single fit step for the specified component.
    
        :param component_name: The name of the component to be fitted.
        """
        param_names = [param for param, (value, fit) in getattr(state,component_name).values_dict.items() if fit]
        def objective_function(params,state,component):
            # Update the parameters in the MolecularSystem
            self.update_parameters(state, params, component)
            # Run DUO and generate output
            self.run_duo()
            # Calculate the RMS error from DUO output
            rms_error = self.output.rms  # Assuming this is correctly calculated and stored
            
            return rms_error
        
        def callback(params):
            rms_error = self.output.rms  # Assuming this is updated each time run_duo() is called
            param_output = ", ".join(f"{name}: {param:.4f}" for name, param in zip(param_names, params))
            #print(f"Current parameters: {param_output} | Current RMS error: {rms_error:.4f}")
            # filename = f'./iteration.txt'
            # self.output.output.to_csv(filename, sep='\t', index=False)
        # Perform the optimization
        
        # Extract parameters for the specified component
        params = self.extract_parameters(state, component_name)
        # Perform optimization using the extracted parameters
        result = minimize(objective_function, params, args=(state, component_name,), callback=callback, method='Nelder-Mead')
        # Update the component's parameters with the optimized values
        self.update_parameters(state, result.x, component_name)
    
    def single_fit_step_off_diag(self, component, verbose=False):
        """
        Performs a single fit step for the specified component.
    
        :param component_name: The name of the component to be fitted.
        """
        param_names = [param for param, (value, fit) in component.values_dict.items() if fit]
        def objective_function(params, component):
            # Update the parameters in the MolecularSystem
            self.update_off_diag(params, component)
            # Run DUO and generate output
            self.run_duo()
            # Calculate the RMS error from DUO output
            rms_error = self.output.rms  # Assuming this is correctly calculated and stored
            
            return rms_error
        
        def callback(params):
            rms_error = self.output.rms  # Assuming this is updated each time run_duo() is called
            param_output = ", ".join(f"{name}: {param:.4f}" for name, param in zip(param_names, params))
            #print(f"Current parameters: {param_output} | Current RMS error: {rms_error:.4f}")
            # filename = f'./iteration.txt'
            # self.output.output.to_csv(filename, sep='\t', index=False)
        # Perform the optimization
        
        # Extract parameters for the specified component
        params = self.extract_off_diag_parameters(component)
        # Perform optimization using the extracted parameters
        result = minimize(objective_function, params, args=(component,), callback=callback, method='Nelder-Mead')
        # Update the component's parameters with the optimized values
        self.update_off_diag(result.x, component)
    
    def fit(self, verbose = False):
        self.fitting_active = True
        # Extract parameters to be varied and their initial values
        initial_params = self.extract_parameters()
        param_names = [param for param, (value, fit) in self.states[0].potential_energy.values_dict.items() if fit]
        # Define the objective function to minimize
        def objective_function(params):
            # Update the parameters in the MolecularSystem
            self.update_parameters(params)
            # Run DUO and generate output
            self.run_duo()
            # Calculate the RMS error from DUO output
            rms_error = self.output.rms  # Assuming this is correctly calculated and stored
            
            return rms_error
        
        def callback(params):
            rms_error = self.output.rms  # Assuming this is updated each time run_duo() is called
            param_output = ", ".join(f"{name}: {param:.4f}" for name, param in zip(param_names, params))
            
            print(f"Current parameters: {param_output} | Current RMS error: {rms_error:.4f}")
            #filename = f'./iteration.txt'
            #self.output.output.to_csv(filename, sep='\t', index=False)
        # Perform the optimization
        result = minimize(objective_function, initial_params, method='Nelder-Mead', callback=callback)


        # Update the system with the optimized parameters
        self.update_parameters(result.x)
        
        if verbose:
            print("Optimization result:")
            print(result.x)
            # Save the DataFrame to a text file        
        # Return the optimization result
        return result
    
    def iterative_fit(self, state_vary=None, off_diagonal_vary=None, convergence_threshold=1e-6, max_iterations=100,verbose=False):
        """
        Iteratively fits the parameters in different blocks until convergence.

        :param params_to_vary: Dictionary of parameters to vary for each component.
        :param convergence_threshold: The RMS change threshold for convergence.
        :param max_iterations: Maximum number of iterations to prevent infinite loops.
        """
        self.fitting_active = True
        previous_rms = float('inf')  # Initialize with a very large number
        for iteration in range(max_iterations):
            # Alternate fitting steps between blocks
            if state_vary:
                for state in state_vary:
                    print(f'Fitting {state}')
                    for component_name in state_vary[state].keys():
                        print(f"FITTING COMPONENT {component_name}")
                        self.single_fit_step(self.states[state], component_name)
                        # After fitting, run DUO to get the new RMS
                        self.run_duo()
                        current_rms = self.output.rms
        
                        # Check convergence criterion
                    rms_change = abs(previous_rms - current_rms)
            if off_diagonal_vary:
                for off_diagonal_name in off_diagonal_vary.keys():
                    print(f"FITTING OFF-DIAGONAL {off_diagonal_name}")
                    self.single_fit_step_off_diag(self.off_diagonal[off_diagonal_name])
                    self.run_duo()
                    current_rms=self.output.rms
            if rms_change < convergence_threshold:
                print(f"Converged after {iteration+1} iterations.")
                print(f"Iteration {iteration+1}: RMS = {current_rms}")
                return
            previous_rms = current_rms
            print(f"Iteration {iteration+1}: RMS = {current_rms}")
            if not state_vary and not off_diagonal_vary:
                raise ValueError("No states or off-diagonals set to vary")

        print("Maximum iterations reached without convergence.")
        if verbose:
            print("Optimization result:")
            print(result.x)        
        # Return the optimization result
        return result