#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  2 20:27:20 2023

@author: shaun
"""
import numpy as np
import pandas as pd
import logging
import matplotlib.pyplot as plt
import matplotlib as mpl
from watchdog.observers import Observer
from watchdog.events import FileSystemEventHandler
import time

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
    def __init__(self,duo_output, output_type='Energies',weights=None,state_number=None):
        self.output_file = duo_output
        if output_type.lower() == 'energies':
            columns = ['##', 'N', 'J', 'p', 'obs', 'calc', 'obs-calc','Weight', '(','States1','vib1','lambda1','sigma1','omega1',')(','states2','vib2', 'lambda2', 'sigma2', 'omega2', ')']
        if output_type.lower() == 'frequencies':
            columns = ['##', 'N', 'J', 'p', 'obs', 'calc', 'obs-calc','Weight', '(','States1','vib1','lambda1','sigma1','omega1',')(','states2','vib2', 'lambda2', 'sigma2', 'omega2', ')','name']
        with open(duo_output) as output_file:
            self.df = pd.read_table(duo_output, delim_whitespace = True, skiprows=6, names=columns)
            self.output = self.df[self.df['obs']!=0.0000].reset_index(drop=True)
            try:
                if state_number:
                    self.rms = self.weighted_rms_state(state_number)
                else:
                    self.rms = self.weighted_rms()
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
        obs = self.output['obs']
        calc = self.output['calc']
        squared_errors = (obs - calc) ** 2
        weighted_squared_errors = squared_errors * self.output['Weight']
        weighted_mse = np.sum(weighted_squared_errors) / np.sum(self.output['Weight'])
        rms = np.sqrt(weighted_mse)
        return rms
    
    def weighted_rms_state(self, state_number):
        """
        Calculate the weighted root mean square error between observed and calculated values.
        """
        filtered_df = self.output[self.output['States1']==state_number]
        obs = filtered_df['obs']
        calc = filtered_df['calc']
        squared_errors = (obs - calc) ** 2
        weighted_squared_errors = squared_errors * filtered_df['Weight']
        weighted_mse = np.sum(weighted_squared_errors) / np.sum(filtered_df['Weight'])
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
            current_block = {'A_values': [],
                             'B_values': []}
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
            elif key.startswith('B') and key != 'BINF':
                current_block['B_values'].append(value)
            else:
                current_block[key] = value

    if current_block:
        block_name = current_block.get('name', '').strip('"')
        blocks.setdefault(current_block_type, {})[block_name] = current_block

    return blocks

def plot_residuals(output_filepath=None,rot_state=False):
    mpl.rcParams['font.size']=10
    plt.rcParams['axes.facecolor'] = 'white'
    data_list_with_J = []
    # Filter out rows where Obs. is 0.0000
    with open(output_filepath, 'r') as f:
        for line in f:
            # Skip header lines or empty lines
            if line.startswith('|') or line.strip() == '':
                continue
            
            # Split the line into columns based on whitespace
            columns = line.split()
            
            # Extract relevant columns: 'J', 'Obs.', 'Calc.', 'v', 'p'
            try:
                J = float(columns[2])  # The 'J' values are actually in the third column
                obs = float(columns[4])
                calc = float(columns[5])
                obs_calc = float(columns[6])
                state = int(columns[9])
                vib = int(columns[10])
                p = columns[3][-1]  # Extract the last character from the 'J p' column
                
                # Include only rows where Obs. is not 0.0000
                if obs != 0.0000:
                    data_list_with_J.append({'J': J, 'Obs.': obs, 'Calc': calc, 'obs-calc': obs_calc, 'state': state, 'vib': vib, 'p': p})
                    
            except (IndexError, ValueError) as e:
                # Handle lines that are not in the expected format
                continue

    # Convert the list of dictionaries to a DataFrame
    df_new_filtered_with_J = pd.DataFrame(data_list_with_J)

    # Determine the number of unique states
    num_states = df_new_filtered_with_J['state'].nunique()

    # Create a figure with a subplot for each state
    fig, axes = plt.subplots(nrows=num_states, ncols=1, figsize=(12, 6 * num_states))


    # Check if there is only one state (axes will not be an array in this case)
    if num_states == 1:
        axes = [axes]

    # Loop through each 'state' group
    for state_idx, (state, state_group) in enumerate(df_new_filtered_with_J.groupby('state')):
        ax = axes[state_idx]  # Get the corresponding subplot axis
        vib_colors = plt.cm.Paired(np.linspace(0, 1, state_group['vib'].nunique()))  # Color map for vibes within a state

        # Loop through each 'vib' group within the current 'state'
        for vib_idx, (vib, group_data) in enumerate(state_group.groupby('vib')):
            # Separate the data by 'p' value
            group_data_p_plus = group_data[group_data['p'] == '+']
            group_data_p_minus = group_data[group_data['p'] == '-']
            jjpp1 = group_data_p_plus['J']*(group_data_p_plus['J']+1)
            jjmp1 = group_data_p_minus['J']*(group_data_p_minus['J']+1)
            # Plot the data with different markers for '+' and '-'
            if rot_state:
                if state==1:
                    ax.scatter(jjpp1, group_data_p_plus['Calc']-3.6334432*jjpp1, color=vib_colors[vib_idx], marker='+', label=f'vib={vib}, p=+', s=90)
                    ax.scatter(jjmp1, group_data_p_minus['Calc']-3.6334432*jjmp1, color=vib_colors[vib_idx], marker='_', label=f'vib={vib}, p=-', s=90)
                elif state==2:
                    ax.scatter(jjpp1, group_data_p_plus['Calc']-3.805503*jjpp1, color=vib_colors[vib_idx], marker='+', label=f'vib={vib}, p=+', s=90)
                    ax.scatter(jjmp1, group_data_p_minus['Calc']-3.805503*jjmp1, color=vib_colors[vib_idx], marker='_', label=f'vib={vib}, p=-', s=90)
                elif state==3:
                    ax.scatter(jjpp1, group_data_p_plus['Calc']-3.6*jjpp1, color=vib_colors[vib_idx], marker='+', label=f'vib={vib}, p=+', s=90)
                    ax.scatter(jjmp1, group_data_p_minus['Calc']-3.6*jjmp1, color=vib_colors[vib_idx], marker='_', label=f'vib={vib}, p=-', s=90)
                elif state==4:
                    ax.scatter(jjpp1, group_data_p_plus['Calc']-3.6*jjpp1, color=vib_colors[vib_idx], marker='+', label=f'vib={vib}, p=+', s=90)
                    ax.scatter(jjmp1, group_data_p_minus['Calc']-3.6*jjmp1, color=vib_colors[vib_idx], marker='_', label=f'vib={vib}, p=-', s=90)
            else:
                ax.scatter(group_data_p_plus['J'], group_data_p_plus['obs-calc'], color=vib_colors[vib_idx], marker='+', label=f'vib={vib}, p=+', s=90)
                ax.scatter(group_data_p_minus['J'], group_data_p_minus['obs-calc'], color=vib_colors[vib_idx], marker='_', label=f'vib={vib}, p=-', s=90)

        ax.set_xlabel('J')
        ax.set_ylabel('Obs-Calc.')
        ax.set_title(f'State {state} obs-calc')
        ax.legend(fontsize=8)
        ax.grid(True)

    # Adjust layout
    plt.tight_layout()
    plt.show()
    
def plot_calc(output_filepath=None,rot_state=False):
    mpl.rcParams['font.size']=10
    plt.rcParams['axes.facecolor'] = 'white'
    data_list_with_J = []
    # Filter out rows where Obs. is 0.0000
    with open(output_filepath, 'r') as f:
        for line in f:
            # Skip header lines or empty lines
            if line.startswith('|') or line.strip() == '':
                continue
            
            # Split the line into columns based on whitespace
            columns = line.split()
            
            # Extract relevant columns: 'J', 'Obs.', 'Calc.', 'v', 'p'
            try:
                J = float(columns[2])  # The 'J' values are actually in the third column
                obs = float(columns[4])
                calc = float(columns[5])
                obs_calc = float(columns[6])
                state = int(columns[9])
                vib = int(columns[10])
                p = columns[3][-1]  # Extract the last character from the 'J p' column
                
                # Include only rows where Obs. is not 0.0000
                if vib<3:
                    data_list_with_J.append({'J': J, 'Obs.': obs, 'Calc': calc, 'obs-calc': obs_calc, 'state': state, 'vib': vib, 'p': p})
                    
            except (IndexError, ValueError) as e:
                # Handle lines that are not in the expected format
                continue

    # Convert the list of dictionaries to a DataFrame
    df_new_filtered_with_J = pd.DataFrame(data_list_with_J)

    # Determine the number of unique states
    num_states = df_new_filtered_with_J['state'].nunique()

    # Create a figure with a subplot for each state
    fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(12, 6 * num_states))


    # Check if there is only one state (axes will not be an array in this case)
    if num_states == 1:
        axes = [axes]

    # Loop through each 'state' group
    for state_idx, (state, state_group) in enumerate(df_new_filtered_with_J.groupby('state')):
        ax = axes  # Get the corresponding subplot axis
        vib_colors = plt.cm.Paired(np.linspace(0, 1, state_group['vib'].nunique()))  # Color map for vibes within a state

        # Loop through each 'vib' group within the current 'state'
        for vib_idx, (vib, group_data) in enumerate(state_group.groupby('vib')):
            # Separate the data by 'p' value
            group_data_p_plus = group_data[group_data['p'] == '+']
            group_data_p_minus = group_data[group_data['p'] == '-']
            jjpp1 = group_data_p_plus['J']*(group_data_p_plus['J']+1)
            jjmp1 = group_data_p_minus['J']*(group_data_p_minus['J']+1)
            # Plot the data with different markers for '+' and '-'
            if rot_state:
                if state==1:
                    ax.scatter(jjpp1, group_data_p_plus['Calc']-3.6334432*jjpp1, color=vib_colors[vib_idx], marker='+', label=f'vib={vib}, p=+', s=20)
                    ax.scatter(jjmp1, group_data_p_minus['Calc']-3.6334432*jjmp1, color=vib_colors[vib_idx], marker='_', label=f'vib={vib}, p=-', s=20)
                elif state==2:
                    ax.scatter(jjpp1, group_data_p_plus['Calc']-3.805503*jjpp1, color=vib_colors[vib_idx], marker='+', label=f'vib={vib}, p=+', s=20)
                    ax.scatter(jjmp1, group_data_p_minus['Calc']-3.805503*jjmp1, color=vib_colors[vib_idx], marker='_', label=f'vib={vib}, p=-', s=20)
                elif state==3:
                    ax.scatter(jjpp1, group_data_p_plus['Calc']-3.6*jjpp1, color=vib_colors[vib_idx], marker='+', label=f'vib={vib}, p=+', s=20)
                    ax.scatter(jjmp1, group_data_p_minus['Calc']-3.6*jjmp1, color=vib_colors[vib_idx], marker='_', label=f'vib={vib}, p=-', s=20)
                elif state==4:
                    ax.scatter(jjpp1, group_data_p_plus['Calc']-3.6*jjpp1, color=vib_colors[vib_idx], marker='+', label=f'vib={vib}, p=+', s=20)
                    ax.scatter(jjmp1, group_data_p_minus['Calc']-3.6*jjmp1, color=vib_colors[vib_idx], marker='_', label=f'vib={vib}, p=-', s=20)
            else:
                ax.scatter(group_data_p_plus['J'], group_data_p_plus['Calc'], color=vib_colors[vib_idx], marker='+', label=f'vib={vib}, p=+', s=20)
                ax.scatter(group_data_p_minus['J'], group_data_p_minus['Calc'], color=vib_colors[vib_idx], marker='_', label=f'vib={vib}, p=-', s=20)

        ax.set_xlabel('J')
        ax.set_ylabel('Calc.')
        ax.set_title(f'State {state} obs-calc')
        ax.legend(fontsize=8)
        ax.grid(True)

    # Adjust layout
    plt.tight_layout()
    plt.show()    

def plot_residuals_watcher(output_filepath=None,rot_state=False):
    mpl.rcParams['font.size']=10
    plt.rcParams['axes.facecolor'] = 'white'
    data_list_with_J = []
    # Filter out rows where Obs. is 0.0000
    with open(output_filepath, 'r') as f:
        for line in f:
            # Skip header lines or empty lines
            if line.startswith('|') or line.strip() == '':
                continue
            
            # Split the line into columns based on whitespace
            columns = line.split()
            
            # Extract relevant columns: 'J', 'Obs.', 'Calc.', 'v', 'p'
            try:
                J = float(columns[2])  # The 'J' values are actually in the third column
                obs = float(columns[4])
                calc = float(columns[5])
                obs_calc = float(columns[6])
                state = int(columns[9])
                vib = int(columns[10])
                p = columns[3][-1]  # Extract the last character from the 'J p' column
                
                # Include only rows where Obs. is not 0.0000
                if obs != 0.0000:
                    data_list_with_J.append({'J': J, 'Obs.': obs, 'Calc': calc, 'obs-calc': obs_calc, 'state': state, 'vib': vib, 'p': p})
                    
            except (IndexError, ValueError) as e:
                # Handle lines that are not in the expected format
                continue
            
    # Convert the list of dictionaries to a DataFrame
    df_new_filtered_with_J = pd.DataFrame(data_list_with_J)

    # Determine the number of unique states
    num_states = df_new_filtered_with_J['state'].nunique()
    # Create a figure with a subplot for each state
    fig, axes = plt.subplots(nrows=num_states, ncols=1, figsize=(12, 6 * num_states))
    if num_states == 1:
        axes = [axes]
    
    def update_plot(df_new_filtered_with_J, fig, axes):
        # Check if there is only one state (axes will not be an array in this case)

        # Loop through each 'state' group
        for state_idx, (state, state_group) in enumerate(df_new_filtered_with_J.groupby('state')):
            ax = axes[state_idx]  # Get the corresponding subplot axis
            vib_colors = plt.cm.Paired(np.linspace(0, 1, state_group['vib'].nunique()))  # Color map for vibes within a state

            # Loop through each 'vib' group within the current 'state'
            for vib_idx, (vib, group_data) in enumerate(state_group.groupby('vib')):
                # Separate the data by 'p' value
                group_data_p_plus = group_data[group_data['p'] == '+']
                group_data_p_minus = group_data[group_data['p'] == '-']
                jjpp1 = group_data_p_plus['J']*(group_data_p_plus['J']+1)
                jjmp1 = group_data_p_minus['J']*(group_data_p_minus['J']+1)
                # Plot the data with different markers for '+' and '-'
                if rot_state:
                    if state==1:
                        ax.scatter(group_data_p_plus['J'], group_data_p_plus['Calc']-3.6334432*jjpp1, color=vib_colors[vib_idx], marker='+', label=f'vib={vib}, p=+', s=20)
                        ax.scatter(group_data_p_minus['J'], group_data_p_minus['Calc']-3.6334432*jjmp1, color=vib_colors[vib_idx], marker='_', label=f'vib={vib}, p=-', s=20)
                    elif state==2:
                        ax.scatter(group_data_p_plus['J'], group_data_p_plus['Calc']-3.805503*jjpp1, color=vib_colors[vib_idx], marker='+', label=f'vib={vib}, p=+', s=20)
                        ax.scatter(group_data_p_minus['J'], group_data_p_minus['Calc']-3.805503*jjmp1, color=vib_colors[vib_idx], marker='_', label=f'vib={vib}, p=-', s=20)
                    elif state==3:
                        ax.scatter(group_data_p_plus['J'], group_data_p_plus['Calc']-3.805503*jjpp1, color=vib_colors[vib_idx], marker='+', label=f'vib={vib}, p=+', s=20)
                        ax.scatter(group_data_p_minus['J'], group_data_p_minus['Calc']-3.805503*jjmp1, color=vib_colors[vib_idx], marker='_', label=f'vib={vib}, p=-', s=20)
                else:
                    ax.scatter(group_data_p_plus['J'], group_data_p_plus['obs-calc'], color=vib_colors[vib_idx], marker='+', label=f'vib={vib}, p=+', s=20)
                    ax.scatter(group_data_p_minus['J'], group_data_p_minus['obs-calc'], color=vib_colors[vib_idx], marker='_', label=f'vib={vib}, p=-', s=20)

            ax.set_xlabel('J')
            ax.set_ylabel('Obs-Calc.')
            ax.set_title(f'State {state} obs-calc')
            ax.legend(fontsize=8)
            ax.grid(True)
        for ax in axes:
            ax.clear()
        fig.canvas.draw()
        fig.canvas.flush_events()
    
    
    class FileChangeHandler(FileSystemEventHandler):
        def __init__(self, df, filename, func):
            self.filename = filename
            self.func = func
            self.df = df
        def on_modified(self, event):
            if event.src_path == self.filename:
                try:
                    self.func(self.df, self.fig, self.axes)
                except:
                    pass
    
    file_path = output_filepath
    
    event_handler = FileChangeHandler(df_new_filtered_with_J, file_path, update_plot(df_new_filtered_with_J, fig, axes))
    observer = Observer()
    observer.schedule(event_handler, path=file_path, recursive=False)
    observer.start()    

    # Keep the script running
    try:
        while True:
            time.sleep(1)
    except KeyboardInterrupt:
        observer.stop()
        print("Stopped file monitoring.")
    
    observer.join()

def pot_file_dataframes(file_path):

# Let's try parsing the actual file while skipping the first three rows
# and then identifying the state name keywords to separate the data.

    # Re-opening the file
    with open(file_path, 'r') as file:
        lines = file.readlines()
    
    # Skipping the first three rows
    lines = lines[3:]
    
    # Initialize the dictionary to hold dataframes for each state
    state_dfs = {}
    current_state = None
    temp_data = []
    
    # Columns as per the user's request
    columns = ['state_number', 'r', 'abinitio', 'calc', 'abinitio-calc', 'weight']
    
    for line in lines:
        # Check if the line is a state keyword
        if 'state' not in line.lower() and line.strip() and not line.strip()[0].isdigit():  # State identifier found
            # If there's previous data, convert it to a DataFrame
            if temp_data:
                state_dfs[current_state] = pd.DataFrame(temp_data, columns=columns)
                temp_data = []  # Reset temp data for next state
            # Update the current state
            current_state = line.strip().split()[0]
        else:
            # Data line, split and append to temp data
            values = line.strip().split(maxsplit=len(columns)-1)
            if values and values[0].isdigit():  # Ensure it's a data line
                temp_data.append(values)
    
    # Convert the last state's data to a DataFrame
    if temp_data:
        state_dfs[current_state] = pd.DataFrame(temp_data, columns=columns)
    
    # Displaying the first few rows of each state DataFrame as an example
    return state_dfs

def plot_abinitio_calc(data_dict):
    """
    Plot abinitio and calc against r for each state in separate axes.
    Skips any states where the character '|' is in the name.

    :param data_dict: Dictionary containing DataFrames for each state.
    """
    # Filter states to exclude those with '|' in their name
    filtered_states = {state: df for state, df in data_dict.items()}#if '|' not in state}

    # Number of states to plot
    num_states = len(filtered_states)
    
    # Create subplots
    fig, axes = plt.subplots(num_states, 1, figsize=(10, 5 * num_states))
    
    # Check if there is only one subplot (in case of a single state)
    if num_states == 1:
        axes = [axes]

    # Plot data for each state
    for ax, (state, df) in zip(axes, filtered_states.items()):
        # Convert data to numeric, as it might be read as strings
        df[['r', 'abinitio', 'calc']] = df[['r', 'abinitio', 'calc']].apply(pd.to_numeric, errors='coerce')
        
        # Plotting
        ax.plot(df['r'], df['abinitio'], label='ab initio', marker='o')
        ax.plot(df['r'], df['calc'], label='calc', marker='x')
        
        # Set the title and labels
        ax.set_title(f'State: {state}')
        ax.set_xlabel('r')
        ax.set_ylabel('Energy')
        ax.legend()

    plt.tight_layout()
    return plt
