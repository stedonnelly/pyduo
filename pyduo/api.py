#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  2 20:27:20 2023

@author: shaun
"""
import numpy as np
import pandas as pd
import logging

import colorlog

pyduo_logger = logging.getLogger('pyduo')

def set_logger(level: int = logging.INFO, filename: str = None):
    """Initialises Python logging, formatting it nicely,
    and optionally printing to a file.
    """
    log_format = '%(asctime)s - ' '%(funcName)s - ' '%(levelname)s - ' '%(message)s'
    bold_seq = '\033[1m'
    colorlog_format = f'{bold_seq} ' '%(log_color)s ' f'{log_format}'
    colorlog.basicConfig(format=colorlog_format)
    pyduo_logger.setLevel(level)

    if filename is not None:
        fh = logging.FileHandler(filename)
        fh.setLevel(level)
        formatter = logging.Formatter(log_format)
        fh.setFormatter(formatter)
        pyduo_logger.addHandler(fh)

class duo_output_object:
    def __init__(self,duo_output, output_type='Energies',weights=None):
        if output_type.lower() == 'energies':
            columns = ['##', 'N', 'J', 'p', 'obs', 'calc', 'obs-calc',' weight', '(','States1','vib1','lambda1','sigma1','omega1',')(','states2','vib2', 'lambda2', 'sigma2', 'omega2', ')']
        if output_type.lower() == 'frequencies':
            columns = ['##', 'N', 'J', 'p', 'obs', 'calc', 'obs-calc',' weight', '(','States1','vib1','lambda1','sigma1','omega1',')(','states2','vib2', 'lambda2', 'sigma2', 'omega2', ')','name']
        with open(duo_output) as output_file:
            self.df = pd.read_table(duo_output, delim_whitespace = True, skiprows=6, names=columns)
            self.output = self.df[self.df['obs']!=0.0000].reset_index(drop=True)
            try:
                if weights.any:
                    self.merged_df = pd.merge(self.output, weights, left_on='obs', right_on='Energy')
                    self.rms = self.weighted_rms()               
            except:
                try:
                    self.rms = self.calculate_rms()
                except TypeError:
                    self.rms = np.nan

    def calculate_rms(self):
        """Calculates the RMS of the 'obs-calc' column."""
        rms = np.sqrt(np.mean(np.square(self.output['obs-calc'])))
        return rms
    
    def weighted_rms(self):
        """
        Calculate the weighted root mean square error between observed and calculated values.
        """
        #self.weight_df = 
        obs = self.merged_df['obs']
        calc = self.merged_df['calc']
        squared_errors = (obs - calc) ** 2
        weighted_squared_errors = squared_errors * self.merged_df['Weight']
        weighted_mse = np.sum(weighted_squared_errors) / np.sum(self.merged_df['Weight'])
        rms = np.sqrt(weighted_mse)
        return rms
    
def recover_all_params(states):
    output = {}
    components = ['potential_energy', 'spin_rot', 'bob_rot', 'lambda_q', 'lambda_p2q','spin_orbit_coupling']
    for key, state in states.items():
        output[key] = {}
        for component in components:
            try:
                output[key][component] = getattr(state, component).output_params_dict()
            except AttributeError:
                pass
    return output

def process_blocks(input_file):
    with open(input_file, 'r') as file:
        file_lines = file.readlines()
    block_keywords = ['poten', 'spin-orbit', 'spin-orbit-x', 'bob-rot', 'spin-rot', 'lambda-q', 'lambda-p2q']
    end_keyword = 'end'
    blocks = {}
    current_block = None
    current_block_type = ""

    for line in file_lines:
        words = line.strip().split()

        if not words:
            continue

        if words[0].lower() in block_keywords:
            if current_block:
                # Use the 'name' value as the key for the block
                block_name = current_block.get('name', '').strip('"')
                blocks.setdefault(current_block_type, {})[block_name] = current_block
            current_block_type = words[0].lower()
            current_block = {'A_values': []}
        elif words[0].lower() == end_keyword:
            if current_block:
                block_name = current_block.get('name', '').strip('"')
                blocks.setdefault(current_block_type, {})[block_name] = current_block
                current_block = None
        elif current_block is not None:
            key = words[0]
            value = words[1] if len(words) > 1 else None
            try:
                value = float(value)
            except (ValueError, TypeError):
                pass
            if key.startswith('A') and key != 'AINF':
                current_block['A_values'].append(value)
            else:
                current_block[key] = value

    if current_block:
        block_name = current_block.get('name', '').strip('"')
        blocks.setdefault(current_block_type, {})[block_name] = current_block

    return blocks
