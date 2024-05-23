#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  2 13:50:53 2023

@author: shaun
"""

import pandas as pd

class EMO:
    def __init__(self, **kwargs):
        # Here we can define attributes specific to EMO functional
        self.V0 = kwargs.get('V0', 0.0)
        self.RE = kwargs.get('RE', 0.0)
        self.DE = kwargs.get('DE', 0.0)
        self.RREF = kwargs.get('RREF', -1.0)
        self.PL = kwargs.get('PL', 0.0)
        self.PR = kwargs.get('PR', 0.0)
        self.NL = kwargs.get('NL', 0.0)
        self.NR = kwargs.get('NR', 0.0)
        self.A_values = kwargs.get('A_values', [])

    def to_dict(self):
        """Converts the EMO attributes to a dictionary"""
        values_dict = {
            "V0": (self.V0, False),
            "RE": (self.RE, False),
            "DE": (self.DE, False),
            "RREF": (self.RREF, False),
            "PL": (self.PL, False),
            "PR": (self.PR, False),
            "NL": (self.NL, False),
            "NR": (len(self.A_values)-1, False),
        }
        for i, val in enumerate(self.A_values):
            values_dict[f"A{i}"] = (val, False)
        return values_dict
    
class BOBLEROY:
    def __init__(self, **kwargs):
        # Here we can define attributes specific to BOBLEROY functional
        self.RE = kwargs.get('RE', 0.0)
        self.RREF = kwargs.get('RREF', -1.0)
        self.P = kwargs.get('P', 0.0)
        self.A_values = kwargs.get('A_values', [])
        self.AINF = kwargs.get('AINF', 0.0)

    def to_dict(self):
        """Converts the BOBLEROY attributes to a dictionary with each value as a tuple
        where the first element is the value and the second element is a boolean indicating
        if the parameter should be varied during fitting."""
        values_dict = {
            "RE": (self.RE, False),
            "RREF": (self.RREF, False),
            "P": (self.P, False),
            "NT": (len(self.A_values)-1, False),  # NT is the count of A_values
        }
        for i, val in enumerate(self.A_values):
            values_dict[f"A{i}"] = (val, False)
        values_dict['AINF'] = (self.AINF, False)
        return values_dict
    
class BOBLEROY_DAMP:
    def __init__(self, **kwargs):
        self.RE = kwargs.get('RE', 0.0)
        self.RREF = kwargs.get('RREF', -1.0)
        self.P = kwargs.get('P', 0.0)
        self.A_values = kwargs.get('A_values', [])
        self.AINF = kwargs.get('AINF', 0.0)
        self.tdamp = kwargs.get('tdamp',0.0)
        self.r0 = kwargs.get('r0',0.0)
        self.alpha = kwargs.get('alpha',0.0)
    
    def to_dict(self):
        """Converts the BOBLEROY attributes to a dictionary with each value as a tuple
        where the first element is the value and the second element is a boolean indicating
        if the parameter should be varied during fitting."""
        values_dict = {
            "RE": (self.RE, False),
            "RREF": (self.RREF, False),
            "P": (self.P, False),
            "NT": (len(self.A_values)-1, False),  # NT is the count of A_values
        }
        for i, val in enumerate(self.A_values):
            values_dict[f"A{i}"] = (val, False)
        values_dict['AINF'] = (self.AINF, False)
        values_dict['tdamp'] = (self.tdamp, False)
        values_dict['r0'] = (self.r0, False)
        values_dict['alpha'] = (self.alpha, False)
        return values_dict
    
class BOBNA:
    def __init__(self, **kwargs):
        self.RE = kwargs.get('RE',0.0)
        self.Maref = kwargs.get('Maref', 0.0)
        self.Ma = kwargs.get('Ma', 0.0)
        self.Mbref = kwargs.get('Mbref',0.0)
        self.Mb = kwargs.get('Mb',0.0)
        self.P = kwargs.get('P',0.0)
        self.ainf = kwargs.get('AINF',0.0)
        self.binf = kwargs.get('BINF',0.0)
        self.A_values = kwargs.get('A_values', [])
        self.B_values = kwargs.get('B_values', [])
        
    def to_dict(self):
        """Converts the BOBNA attributes to a dictionary with each value as a tuple
        where the first element is the value and the second element is a boolean indicating
        if the parameter should be varied during fitting."""
        
        values_dict = {
            "RE": (self.RE, False),
            "P": (self.P, False),
            "NTa": (len(self.A_values)-1, False),  # NTa is the count of A_values
            "NTb": (len(self.A_values)-1, False),  # NTb is the count of B_values
        }
        
        for i, val in enumerate(self.A_values):
            values_dict[f'A{i}'] = (val, False)
        for i, val in enumerate(self.B_values):
            values_dict[f'B{i}'] = (val, False)
        values_dict['Maref'] = (self.Maref, False)
        values_dict['Ma'] = (self.Ma, False)
        values_dict['Mbref'] = (self.Mbref, False)
        values_dict['Mb'] = (self.Mb, False)
        values_dict['AINF'] = (self.ainf, False)
        values_dict['BINF'] = (self.binf, False)
        return values_dict
        
class POLYNOM_DECAY_24:
    def __init__(self, **kwargs):
        self.RE = kwargs.get('RE', 0.0)
        self.BETA = kwargs.get('BETA', 0.0)
        self.GAMMA = kwargs.get('GAMMA', 0.0)
        self.P = kwargs.get('P', 0.0)
        self.A_values = kwargs.get('A_values')
        self.AINF = kwargs.get('AINF', 0.0)
        
    def to_dict(self):
        
        values_dict = {
            "RE": (self.RE, False),
            "GAMMA": (self.GAMMA, False),
            "BETA": (self.BETA, False),
            "P": (self.P, False)}
        for i, val in enumerate(self.A_values):
            values_dict[f'A{i}'] = (val, False)
        values_dict['AINF'] = (self.AINF, False)
        return values_dict
        
class LORENTZ:
    def __init__(self, **kwargs):
        self.V0 = kwargs.get('V0',0.0)
        self.RE = kwargs.get('RE',0.0)
        self.GAMMA = kwargs.get('GAMMA',0.0)
        self.A_values = kwargs.get('A_values',0.0)
    
    def to_dict(self):
        
        values_dict = {
            "V0": (self.V0, False),
            "RE": (self.RE, False),
            "GAMMA": (self.GAMMA, False)
            }
        for i, val in enumerate(self.A_values):
            values_dict[f'A{i}'] = (val, False)
        return values_dict
        
class GRID:
    def __init__(self, data_csv):
        self.values_df = pd.read_csv(data_csv, names=['R', 'VALUES'],skiprows=1)
        self.R = self.values_df['R']
        self.VALUES = self.values_df['VALUES']
        
    def to_dict(self):
        values_dict = self.values_df.set_index(self.values_df.columns[0])[self.values_df.columns[1]].to_dict()
        return values_dict
        
class MLR:
    def __init__(self,**kwargs):
        # Here we can define attributes specific to EMO functional
        self.V0 = kwargs.get('V0', 0.0)
        self.RE = kwargs.get('RE', 0.0)
        self.DE = kwargs.get('DE', 0.0)
        self.RREF = kwargs.get('RREF', -1.0)
        self.P = kwargs.get('P', 0.0)
        self.NL = kwargs.get('NL', 0.0)
        self.NR = kwargs.get('NR', 0.0)
        self.A_values = kwargs.get('A_values', [])
        self.B_values = kwargs.get('B_values', [])
        
    def to_dict(self):
        """Converts the EMO attributes to a dictionary"""
        values_dict = {
            "V0": (self.V0, False),
            "RE": (self.RE, False),
            "DE": (self.DE, False),
            "RREF": (self.RREF, False),
            "P": (self.P, False),
            "NL": (self.NL, False),
            "NR": (len(self.B_values)-1, False),
        }
        for i, val in enumerate(self.B_values):
            values_dict[f"B{i}"] = (val, False)
        for i, val in enumerate(self.A_values):
            values_dict[f"A{i+1}"] = (val, False)
        return values_dict


class POLYNOM:
    def __init__(self, **kwargs):
        self.A_values = kwargs.get('A_values', [])
        self.RE = kwargs.get('RE', 0.0)
        
    def to_dict(self):
        values_dict = {}
        
        if self.RE:
            values_dict['RE'] = (self.RE, False)
        for i, val in enumerate(self.A_values):
            values_dict[f"A{i+1}"] = (val, False)
        return values_dict


















