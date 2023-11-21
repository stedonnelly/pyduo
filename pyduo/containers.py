from .functionals import EMO
from mendeleev import element
import csv
import pandas as pd
from .components import *

class MoleculeInformation:
    def __init__(self, atom1, atom2, mass1, mass2):
        self.atom1 = atom1
        self.atom2 = atom2
        self.mass1 = mass1
        self.mass2 = mass2
    
    def __str__(self):
        lines = [
            f"masses {self.mass1:.4f} {self.mass2:.4f}",
            f"molecule {self.atom1}{self.atom2}"
        ]
        return "\n".join(lines)

class Grid:
    def __init__(self, npoints, range_start, range_end, grid_type):
        self.npoints = npoints
        self.range_start = range_start
        self.range_end = range_end
        self.grid_type = grid_type
    
    def __str__(self):
        lines = [
            "grid",
            f"  npoints {self.npoints}",
            f"  range  {self.range_start:.2f}, {self.range_end:.2f}",
            f"  type {self.grid_type}",
            "end"
        ]
        return "\n".join(lines)

class ContractionBlock:
    def __init__(self, potential_blocks=[]):
        """Initializes the ContractionBlock with a list of PotentialEnergyBlock objects."""
        self.potential_blocks = potential_blocks
    
    def extract_DE(self,block):
        try:
            de = block.values_dict['DE'][0]
            return de
        except:
            key = list(block.values_dict.keys())[-1]
            de = block.values_dict[key]
            return de
    
    def __str__(self):
        DE_values = [self.extract_DE(block) for block in self.potential_blocks]
        DE_values_str = " ".join([f"{DE:.14E}" for DE in DE_values])
        return f"CONTRACTION\n    enermax {DE_values_str}\nEND"

class DiagonalizerBlock:
    def __init__(self):
        self.type_of_object = "Diagonalizer"
    def __str__(self):
        return f"DIAGONALIZER\n    SYEV\nend"


class Fitting:
    def __init__(self, jlist, output_name, csv_path):
        self.jlist = jlist
        self.output_name = output_name
        self.csv_path = csv_path
        self.df = pd.read_csv(self.csv_path)
        self.weights = self.df.Weight

    def _read_csv(self):
        """Reads the energies from the CSV file."""
        self.df = pd.read_csv(self.csv_path)

    def __str__(self):
        energies_str = self.df.to_csv(sep='\t', index=False,header=False)
        lines = [
            "FITTING",
            f"JLIST\t{', '.join(self.jlist)}",
            "itmax 0",
            "fit_factor\t1e12",
            "lock\t10000",
            f"output\t{self.output_name}",
            "robust 0.0000001",
            "target_rms 0.000001",
            "energies",
            energies_str[:-1],
            "end"
        ]
        return "\n".join(lines)

class FreqFitting:
    def __init__(self, jlist, output_name, csv_path):
        self.jlist = jlist
        self.output_name = output_name
        self.csv_path = csv_path
        self.df = pd.read_csv(self.csv_path)
        self.weights = self.df.Weight

    def _read_csv(self):
        """Reads the energies from the CSV file."""
        self.df = pd.read_csv(self.csv_path)

    def __str__(self):
        energies_str = self.df.to_csv(sep='\t', index=False,header=False)
        lines = [
            "FITTING",
            f"JLIST\t{', '.join(self.jlist)}",
            "itmax 0",
            "fit_factor\t1e12",
            "lock\t10000",
            f"output\t{self.output_name}",
            "robust 0.0000001",
            "target_rms 0.000001",
            "frequencies",
            energies_str[:-1],
            "end"
        ]
        return "\n".join(lines)

class AbInitioBlock:
    def __init__(self, block_object, units, fit_factor, exclude_attrs=["expansion_type", "values_dict", "symmetry"]):
        self.block_object = block_object
        self.units = units
        self.fit_factor = fit_factor
        self.values_df = None  # Will hold the pandas DataFrame for ab initio values
        self.exclude_attrs = exclude_attrs

    def set_values_from_dataframe(self, df):
        """
        Sets the ab initio values from a pandas DataFrame.
        """
        self.values_df = df

    def _construct_preamble(self):
        """
        Constructs the preamble using the attributes from the passed object,
        excluding those in the exclude_attrs list.
        """
        if self.block_object.object_type == "poten":
            preamble = [f"abinitio poten {getattr(self.block_object, 'state_number')}"]
            
        if self.block_object.object_type == "spin-orbit":
            preamble = [f"abinitio spin-orbit {getattr(self.block_object, 'state_numbers')[0]} {getattr(self.block_object, 'state_numbers')[1]}"]
        attributes = [
            attr for attr in dir(self.block_object)
            if not callable(getattr(self.block_object, attr)) and not attr.startswith("_") and attr not in self.exclude_attrs
        ]
        for attr in attributes:
            value = getattr(self.block_object, attr)
            if attr == "state_number":
                continue 
            if attr == "lambda_":
                attr = "lambda"
            if attr == "ab_initio_data":
                continue
            preamble.append(f"{attr} {value}")
        return preamble

    def _values_to_str(self):
        """
        Converts the values DataFrame to a string representation.
        """
        
        if self.values_df is not None:
            df_str = self.values_df.to_csv(sep='\t', index=False,header=False)
            return df_str
        return ""

    def generate_abinitio_old(self):
        preamble_str = "\n".join(self._construct_preamble())
        values_str = self._values_to_str()
        lines = [
            preamble_str,
            "type grid",  # This remains static
            f"units {self.units}",
            f"fit_factor {self.fit_factor}",
            "values",
            values_str,
            "end"
        ]
        return "\n".join(lines)
    
    def generate_abinitio(self):
        return str(self.block_object.create_ab_initio_block())