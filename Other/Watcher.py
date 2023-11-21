import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
from watchdog.observers import Observer
from watchdog.events import FileSystemEventHandler
import time
import pyduo

mpl.rcParams['font.size'] = 12
plt.rcParams['axes.facecolor'] = 'white'

# Define your data processing and plotting function here
def process_and_plot_data(filename):
    # Add your data processing logic here
    # For example, reading from a CSV file:
    data = pd.read_csv(filename)
    plt.figure()
    plt.plot(data)  # Replace with your actual plotting code
    plt.show()

class FileChangeHandler(FileSystemEventHandler):
    def __init__(self, filename):
        self.filename = filename

    def on_modified(self, event):
        if event.src_path == self.filename:
            process_and_plot_data(self.filename)

if __name__ == "__main__":
    file_path = './your_data_file.csv'  # Replace with your file path

    # Set up the file watcher
    event_handler = FileChangeHandler(file_path)
    observer = Observer()
    observer.schedule(event_handler, path=file_path, recursive=False)
    observer.start()

    print("Monitoring file for changes. Press Ctrl+C to stop.")

    # Keep the script running
    try:
        while True:
            time.sleep(1)
    except KeyboardInterrupt:
        observer.stop()
        print("Stopped file monitoring.")
    
    observer.join()