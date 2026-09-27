# get_mouse_mapping_pct.py
# makes a csv file of mouse reads as a ratio of human reads

import os
import pandas as pd

# Get all files in the current directory
files = os.listdir("../logs")

# Filter files that don't contain "mouse"
non_mouse_files = [file.replace("_log.bowtie2", "") for file in files 
                  if "mouse" not in file.lower() and "_log.bowtie2" in file]

non_mouse_files.sort()
non_mouse_files[1:3]


def extract_mapped_reads(bowtie_file_path):
    """
    Extract the sum of mapped reads from a bowtie2 log file.
    
    Args:
        bowtie_file_path (str): Path to the bowtie2 log file
        
    Returns:
        int: Sum of mapped reads (concordantly aligned exactly 1 time + >1 times)
        float: Overall alignment rate percentage
    """
    try:
        with open(bowtie_file_path, 'r') as file:
            content = file.read()
            
            # Extract reads aligned exactly 1 time
            aligned_once_line = [line for line in content.split('\n') if "aligned concordantly exactly 1 time" in line][0]
            aligned_once = int(aligned_once_line.split()[0])
            
            # Extract reads aligned more than 1 time
            aligned_multiple_line = [line for line in content.split('\n') if "aligned concordantly >1 times" in line][0]
            aligned_multiple = int(aligned_multiple_line.split()[0])
            
            # Calculate sum of mapped reads
            sum_mapped_reads = aligned_once + aligned_multiple
            
            # Extract overall alignment rate
            alignment_rate_line = [line for line in content.split('\n') if "overall alignment rate" in line][0]
            alignment_rate = float(alignment_rate_line.split('%')[0])
            
            return sum_mapped_reads, alignment_rate
            
    except Exception as e:
        print(f"Error processing file {bowtie_file_path}: {str(e)}")
        return 0, 0.0

# Create a list to store the results
results = []

# Process each file
for file in non_mouse_files:
    try:
        human_reads, human_rate = extract_mapped_reads("logs/" + file + "_log.bowtie2")
        mouse_reads, mouse_rate = extract_mapped_reads("logs/" + file + "_mouse_log.bowtie2")
        
        results.append({
            'file': file,
            'human_reads': human_reads,
            'human_alignment_pct': human_rate,
            'mouse_reads': mouse_reads,
            'mouse_alignment_pct': mouse_rate,
            'human_mouse_ratio': human_reads/mouse_reads/100 if mouse_reads > 0 else float('inf')
        })
    except Exception as e:
        print(f"Error processing {file}: {str(e)}")

# Create a pandas DataFrame
df = pd.DataFrame(results)

df.to_csv('mapping_results.csv', index=False)