#FastQC
import os
import subprocess
import shutil
import glob
import pandas as pd

# List of input files
G_files = ["G_1.fq", "G_2.fq"]

# Create output directory
output_dir = "/home/mga/Seq_data1"
os.makedirs(output_dir, exist_ok=True)

# Check if FastQC is installed
if shutil.which("fastqc") is None:
    raise FileNotFoundError("FastQC is not installed or not in PATH.")

# Run FastQC with output directory specified
try:
    subprocess.run(["fastqc", "-o", output_dir] + G_files, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
except subprocess.CalledProcessError as e:
    print(f"FastQC failed with error:\n{e.stderr.decode()}")
    exit(1)

# Verify that reports were created
for file in G_files:
    base_name = os.path.splitext(os.path.basename(file))[0]  # Extract filename without extension
    html_file = os.path.join(output_dir, f"{base_name}_fastqc.html")

    if os.path.exists(html_file):
        print(f"Report generated: {html_file}")
    else:
        print(f"Warning: {html_file} was not created!")

print("FastQC analysis completed. Reports are in", output_dir)



#-------------------------------------------------------------------



# Set working directory
working_dir = "/home/mga/Seq_data2"
trimmomatic_path = "/home/mga/trimmomatic"

# Define input files and output names
input_file1 = os.path.join(working_dir, "G_1.fq")
input_file2 = os.path.join(working_dir, "G_2.fq")
base_output = os.path.join(working_dir, "SunClipped")

# Ensure the working directory exists
os.makedirs(working_dir, exist_ok=True)

# Check if input files exist
for file in [input_file1, input_file2]:
    if not os.path.exists(file):
        raise FileNotFoundError(f"Error: {file} not found!")

# Ensure Trimmomatic is executable
if not shutil.which(trimmomatic_path) and not os.path.isfile(trimmomatic_path):
    raise FileNotFoundError(f"Error: Trimmomatic not found at {trimmomatic_path}")

# Define Trimmomatic parameters
params = [
    "PE", "-threads", "12",
    input_file1, input_file2,
    "-baseout", base_output,
    "ILLUMINACLIP:Adapter.fa:4:30:10",
    "SLIDINGWINDOW:4:30",
    "MINLEN:30"
]

# Run Trimmomatic
try:
    result = subprocess.run(
        ["bash", trimmomatic_path] + params,
        cwd=working_dir,
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )
    print("Trimmomatic completed successfully.")
    print(result.stdout)
except subprocess.CalledProcessError as e:
    print("Error running Trimmomatic:")
    print(f"STDOUT:\n{e.stdout}")
    print(f"STDERR:\n{e.stderr}")
    exit(1)



#------------------------------------------------------------------------



# Define working directory
working_dir = "/home/mga/Seq_data2"
os.makedirs(working_dir, exist_ok=True)

# Find all tDNA files in the working directory
tDNA_files = glob.glob(os.path.join(working_dir, "tDNA*.fa"))

# Check if any tDNA files were found
if not tDNA_files:
    print("No tDNA files found for indexing.")
    exit(1)

# Check if bwa-mem2 is installed
if shutil.which("bwa-mem2") is None:
    raise FileNotFoundError("Error: bwa-mem2 not found in PATH.")

# Iterate through each file and run bwa-mem2 index command
for tDNA in tDNA_files:
    try:
        subprocess.run(["bwa-mem2", "index", tDNA], check=True)
        print(f"Indexing completed for {tDNA}")
    except subprocess.CalledProcessError as e:
        print(f"Error indexing {tDNA}: {e}")




#----------------------------------------------------------------------------




# Define working directory
working_dir = "/home/mga/Seq_data2"
os.makedirs(working_dir, exist_ok=True)

# Automatically detect reference files
tDNA_files = sorted(glob.glob(os.path.join(working_dir, "tDNA*.fa")))

# Check if bwa-mem2 is installed
if shutil.which("bwa-mem2") is None:
    raise FileNotFoundError("Error: bwa-mem2 not found in PATH.")

# Define input read files
read_1 = os.path.join(working_dir, "SunClipped_1P")
read_2 = os.path.join(working_dir, "SunClipped_2P")

# Ensure input read files exist
for read_file in [read_1, read_2]:
    if not os.path.exists(read_file):
        raise FileNotFoundError(f"Error: Input read file {read_file} not found!")

# Ensure there are reference files
if not tDNA_files:
    raise FileNotFoundError("Error: No tDNA reference files found!")

# Iterate through tDNA reference files and run bwa-mem2
for idx, ref_file in enumerate(tDNA_files, start=1):
    output_sam = os.path.join(working_dir, f"SunT{idx}.sam")

    try:
        subprocess.run(
            ["bwa-mem2", "mem", ref_file, read_1, read_2, "-o", output_sam, "-t", "12"],
            check=True
        )
        print(f"Alignment completed: {ref_file} → {output_sam}")
    except subprocess.CalledProcessError as e:
        print(f"Error running bwa-mem2 for {ref_file}: {e}")




#--------------------------------------------------------------------------------




#Add First 10 lines to first10lines_all.sam
# Define directories
input_dir = "/home/mga/Seq_data2"
output_dir = "/home/mga/Seq_data1"
output_samcheck = os.path.join(output_dir, "first10lines_all.sam")

# Ensure output directory exists
os.makedirs(output_dir, exist_ok=True)

# Detect all SunT*.sam files automatically
sam_files = sorted(glob.glob(os.path.join(input_dir, "SunT*.sam")))

# Check if any SAM files were found
if not sam_files:
    print("No SunT*.sam files found.")
    exit(1)

# Open the output file for writing
with open(output_samcheck, "w") as outfile:
    for sam_file in sam_files:
        try:
            # Read the first 10 lines using Python instead of `head`
            with open(sam_file, "r") as infile:
                for _ in range(10):
                    line = infile.readline()
                    if not line:
                        break  # Stop if file has fewer than 10 lines
                    outfile.write(line)

            print(f"First 10 lines from {sam_file} added to {output_samcheck}")

        except Exception as e:
            print(f"Error processing {sam_file}: {e}")

print(f"All files processed. Output saved in {output_samcheck}.")




#-------------------------------------------------------------------------------




# Define directories
input_dir = "/home/mga/Seq_data2"
output_dir = "/home/mga/Seq_data2"

# Ensure output directory exists
os.makedirs(output_dir, exist_ok=True)

# Automatically detect SunT*.sam files
sam_files = sorted(glob.glob(os.path.join(input_dir, "SunT*.sam")))

# Check if samtools is installed
if shutil.which("samtools") is None:
    raise FileNotFoundError("Error: samtools not found in PATH.")

# Check if there are any SAM files to process
if not sam_files:
    print("No SunT*.sam files found for conversion.")
    exit(1)

# Convert SAM to BAM for each file
for sam_file in sam_files:
    bam_file = os.path.join(output_dir, os.path.basename(sam_file).replace(".sam", ".bam"))
    try:
        # Run the samtools view command to convert SAM to BAM
        subprocess.run(
            ["samtools", "view", "-@12", "-bS", sam_file, "-o", bam_file],
            check=True
        )
        print(f"Conversion successful: {sam_file} -> {bam_file}")
    except subprocess.CalledProcessError as e:
        print(f"Error converting {sam_file} to BAM: {e}")

print(f"All files processed. BAM files saved in {output_dir}.")




#-------------------------------------------------------------------------------------




# Define directories
input_dir = "/home/mga/Seq_data2"
output_dir = "/home/mga/Seq_data2"

# Ensure output directory exists
os.makedirs(output_dir, exist_ok=True)

# Automatically detect SunT*.bam files
bam_files = sorted(glob.glob(os.path.join(input_dir, "SunT*.bam")))

# Check if samtools is installed
if shutil.which("samtools") is None:
    raise FileNotFoundError("Error: samtools not found in PATH.")

# Check if there are any BAM files to process
if not bam_files:
    print("No SunT*.bam files found for sorting.")
    exit(1)

# Sort each BAM file and output sorted BAM
for bam_file in bam_files:
    sorted_bam = os.path.join(output_dir, os.path.basename(bam_file).replace(".bam", "s.bam"))
    
    try:
        # Run the samtools sort command to sort the BAM file
        subprocess.run(
            ["samtools", "sort", "-@12", bam_file, "-o", sorted_bam],
            check=True
        )
        print(f"Sorting successful: {bam_file} -> {sorted_bam}")
    except subprocess.CalledProcessError as e:
        print(f"Error sorting {bam_file}: {e}")

print(f"All BAM files processed. Sorted files saved in {output_dir}.")




#-------------------------------------------------------------------------------------




# Define directories
input_dir = "/home/mga/Seq_data2"
output_dir = "/home/mga/Seq_data2"

# Ensure output directory exists
os.makedirs(output_dir, exist_ok=True)

# Automatically detect SunT*s.bam files
bam_files_to_merge = sorted(glob.glob(os.path.join(input_dir, "SunT*s.bam")))

# Check if samtools is installed
if shutil.which("samtools") is None:
    raise FileNotFoundError("Error: samtools not found in PATH.")

# Ensure there are BAM files to merge
if not bam_files_to_merge:
    print("No SunT*s.bam files found for merging.")
    exit(1)

# Merge BAM files
merged_bam = os.path.join(output_dir, "SunT_merge.bam")
try:
    subprocess.run(
        ["samtools", "merge", merged_bam] + bam_files_to_merge,
        check=True
    )
    print(f"Merge completed successfully: {merged_bam}")
except subprocess.CalledProcessError as e:
    print(f"Error merging BAM files: {e}")

# Remove duplicates
merged_bam_no_duplicates = os.path.join(output_dir, "SunT_merge2.bam")
try:
    subprocess.run(
        ["samtools", "rmdup", "-sS", merged_bam, merged_bam_no_duplicates],
        check=True
    )
    print(f"Duplicate removal completed successfully: {merged_bam_no_duplicates}")
except subprocess.CalledProcessError as e:
    print(f"Error removing duplicates: {e}")

# Sort the BAM file
sorted_bam = os.path.join(output_dir, "SunTs_merge.bam")
try:
    subprocess.run(
        ["samtools", "sort", "-@12", merged_bam_no_duplicates, "-o", sorted_bam],
        check=True
    )
    print(f"Sorting completed successfully: {sorted_bam}")
except subprocess.CalledProcessError as e:
    print(f"Error sorting BAM file: {e}")

print(f"All tasks completed. Final sorted BAM file saved as {sorted_bam}.")





#---------------------------------------------------------------------------------------





# Define directories and filenames
input_dir = "/home/mga/Seq_data2"
output_dir = "/home/mga/Seq_data2"
reference_genome = "Arab_Genome.fa"
read_1 = "SunClipped_1P"
read_2 = "SunClipped_2P"

# Ensure output directory exists
os.makedirs(output_dir, exist_ok=True)

# Check if bwa-mem2 and samtools are installed
if shutil.which("bwa-mem2") is None:
    raise FileNotFoundError("Error: bwa-mem2 not found in PATH.")
if shutil.which("samtools") is None:
    raise FileNotFoundError("Error: samtools not found in PATH.")

# Step 1: Run bwa-mem2 to align the reads to the reference genome
sunG_sam = os.path.join(output_dir, "SunG.sam")
bwa_command = [
    "bwa-mem2", "mem", reference_genome, read_1, read_2,
    "-o", sunG_sam, "-t", "12"
]

try:
    subprocess.run(bwa_command, check=True)
    print(f"Alignment completed successfully: {sunG_sam}")
except subprocess.CalledProcessError as e:
    print(f"Error running bwa-mem2: {e}")

# Step 2: Convert SAM to BAM using samtools view
sunG_bam = os.path.join(output_dir, "SunG.bam")
samtools_view_command = [
    "samtools", "view", "-@12", "-bS", sunG_sam, "-o", sunG_bam
]

try:
    subprocess.run(samtools_view_command, check=True)
    print(f"SAM to BAM conversion completed successfully: {sunG_bam}")
except subprocess.CalledProcessError as e:
    print(f"Error converting SAM to BAM: {e}")

# Step 3: Sort the BAM file using samtools sort
sunGs_bam = os.path.join(output_dir, "SunGs.bam")
samtools_sort_command = [
    "samtools", "sort", "-@12", sunG_bam, "-o", sunGs_bam
]

try:
    subprocess.run(samtools_sort_command, check=True)
    print(f"BAM sorting completed successfully: {sunGs_bam}")
except subprocess.CalledProcessError as e:
    print(f"Error sorting BAM file: {e}")






#-----------------------------------------------------------------------





# Define directories and filenames
input_bam = "/home/mga/Seq_data2/SunTs_merge.bam"
output_dir = "/home/mga/Seq_data2"
output_bam = os.path.join(output_dir, "SunTFlanks.bam")
output_fq = os.path.join(output_dir, "SunTFlanks.fq")
output_sam = os.path.join(output_dir, "SunInsertions.sam")
reference_genome = "Arab_Genome.fa"

# Ensure output directory exists
os.makedirs(output_dir, exist_ok=True)

# Check if samtools and bwa-mem2 are installed
if shutil.which("samtools") is None:
    raise FileNotFoundError("Error: samtools not found in PATH.")
if shutil.which("bwa-mem2") is None:
    raise FileNotFoundError("Error: bwa-mem2 not found in PATH.")

# Step 1: Filter BAM file with samtools view
filter1_command = [
    "samtools", "view", "-@12", "-u", "-f4", "-F8", input_bam, "-o", output_bam
]
try:
    subprocess.run(filter1_command, check=True)
    print(f"BAM filtering completed successfully: {output_bam}")
except subprocess.CalledProcessError as e:
    print(f"Error filtering BAM file: {e}")

# Step 2: Convert filtered BAM file to FASTQ using samtools bam2fq
output_fq = "SunTFlanks.fq"  # Specify the desired output name

newfq_command = ["samtools", "bam2fq", "-@12", "SunTFlanks.bam"]

try:
    with open(output_fq, "w") as fq_file:  # Redirect output to the specified file
        subprocess.run(newfq_command, stdout=fq_file, check=True)
    print(f"BAM to FASTQ conversion completed successfully: {output_fq}")
except subprocess.CalledProcessError as e:
    print(f"Error converting BAM to FASTQ: {e}")


# Step 3: Align reads using bwa-mem2
map2G_command = [
    "bwa-mem2", "mem", reference_genome, output_fq, "-o", output_sam, "-t", "12"
]
try:
    subprocess.run(map2G_command, check=True)
    print(f"Map to Genome alignment completed successfully: {output_sam}")
except subprocess.CalledProcessError as e:
    print(f"Error running Map to Genome: {e}")






#-----------------------------------------------------------------------------------------





# Define directories and filenames
input_sam = "/home/mga/Seq_data2/SunInsertions.sam"
output_dir = "/home/mga/Seq_data2"
output_bam = os.path.join(output_dir, "SunInsertions.bam")
output_sorted_bam = os.path.join(output_dir, "SunInsertions2.bam")
output_csv = os.path.join(output_dir, "SunInsertions.csv")

# Ensure output directory exists
os.makedirs(output_dir, exist_ok=True)

# Check if samtools is installed
if shutil.which("samtools") is None:
    raise FileNotFoundError("Error: samtools not found in PATH.")

# Step 1: Convert SAM to BAM using samtools view
samtools_view_command = [
    "samtools", "view", "-bS", "-@12", input_sam, "-o", output_bam
]
try:
    subprocess.run(samtools_view_command, check=True)
    print(f"SAM to BAM conversion completed successfully: {output_bam}")
except subprocess.CalledProcessError as e:
    print(f"Error converting SAM to BAM: {e}")

# Step 2: Sort the BAM file using samtools sort
samtools_sort_command = [
    "samtools", "sort", "-@12", output_bam, "-o", output_sorted_bam
]
try:
    subprocess.run(samtools_sort_command, check=True)
    print(f"BAM sorting completed successfully: {output_sorted_bam}")
except subprocess.CalledProcessError as e:
    print(f"Error sorting BAM file: {e}")

# Step 3: Convert the sorted BAM file to CSV using samtools view
samtools_view_csv_command = [
    "samtools", "view", "-F4", "-@12", output_sorted_bam
]
try:
    with open(output_csv, "w") as csv_file:
        subprocess.run(samtools_view_csv_command, stdout=csv_file, check=True)
    print(f"BAM to CSV conversion completed successfully: {output_csv}")
except subprocess.CalledProcessError as e:
    print(f"Error converting BAM to CSV: {e}")







#-----------------------------------------------------------------------------------------------------





#Create file called Soft_Clipped_reads.csv and add column 2 and 3 from SunInsertions.csv for rows in which column 5 contains "S".
# Define input and output file paths
input_table = "/home/mga/Seq_data2/SunInsertions.csv"
output_table = "/home/mga/Seq_data2/Soft_clipped_reads.csv"

# Check if the input table exists
if not os.path.exists(input_table):
    raise FileNotFoundError(f"Input table not found: {input_table}")

# Read the table into a pandas DataFrame (assuming CSV format, adjust delimiter if necessary)
df = pd.read_csv(input_table, delimiter='\t', on_bad_lines='skip', header=None)

# Filter rows where column 5 (index 4) contains the letter 'S'
filtered_df = df[df[5].str.contains('S', na=False)]

# Extract columns 2 (index 1) and 3 (index 2)
result_df = filtered_df[[2, 3]]

# Save the filtered result to a new file
result_df.to_csv(output_table, index=False)

print(f"Filtered data saved to {output_table}")

# Define source and destination file paths
source_path = output_table
destination_path = "/home/mga/Seq_data1/Soft_clipped_reads.csv"

# Check if source file exists
if not os.path.exists(source_path):
    raise FileNotFoundError(f"Source file not found: {source_path}")

# Copy the file
shutil.copy(source_path, destination_path)
print(f"File copied from {source_path} to {destination_path}")

# Copy SunInsertions.csv as well
shutil.copy(input_table, "/home/mga/Seq_data1/SunInsertions.csv")
print(f"File copied from {input_table} to /home/mga/Seq_data1/SunInsertions.csv")






#---------------------------------------------------------------------------------------------------------





# index SunGs.bam, view with Soft_Clipped_reads.csv coordinates as (x:y), sort final file
# Define the BAM file you want to index
bam_file = "SunGs.bam"

# Check if BAM file exists
if not os.path.exists(bam_file):
    raise FileNotFoundError(f"BAM file not found: {bam_file}")

# Run the samtools index command
try:
    subprocess.run(["samtools", "index", bam_file], check=True)
    print(f"Indexing of {bam_file} completed successfully.")
except subprocess.CalledProcessError as e:
    print(f"Error indexing BAM file: {e}")

# Define the input CSV file and the output BAM file
SoftClips_csv = "Soft_clipped_reads.csv"
output_Flanks = "Flanks.bam"

# Check if CSV file exists
if not os.path.exists(SoftClips_csv):
    raise FileNotFoundError(f"CSV file not found: {SoftClips_csv}")

# Read the CSV file into a pandas DataFrame (skip first row directly)
df = pd.read_csv(SoftClips_csv, header=None, skiprows=1)  # Skip the first row directly

# Prepare the list of X values in the format column1:column2
X_values = [f"{row[0]}:{row[1]}" for _, row in df.iterrows()]

# Construct the samtools view command
samtools_command = ["samtools", "view", "SunGs.bam"] + X_values + ["-o", output_Flanks]

# Run the samtools command
try:
    subprocess.run(samtools_command, check=True)
    print(f"Samtools command executed: {' '.join(samtools_command)}")
except subprocess.CalledProcessError as e:
    print(f"Error running samtools view: {e}")

# Define input and output BAM file names for sorting
input_Flanks = "Flanks.bam"
output_FlanksS = "FlanksS.bam"

# Run samtools sort
try:
    subprocess.run(["samtools", "sort", input_Flanks, "-o", output_FlanksS], check=True)
    print(f"Sorted BAM file saved as {output_FlanksS}")
except subprocess.CalledProcessError as e:
    print(f"Error sorting BAM file: {e}")




#-------------------------------------------------------------------------------------------------------



#Merge bam files, remove PCR duplicates
# Define file names
output_merged_bam = "MergeInsertions.bam"
sorted_merged_bam = "MergeInsertionsS.bam"
input_mergers = ["FlanksS.bam", "SunInsertions.bam"]

# Check if input files exist
for file in input_mergers:
    if not os.path.exists(file):
        raise FileNotFoundError(f"Input file not found: {file}")

# Run samtools merge
try:
    subprocess.run(["samtools", "merge", "-o", output_merged_bam] + input_mergers, check=True)
    print(f"Merged BAM file saved as {output_merged_bam}")
except subprocess.CalledProcessError as e:
    print(f"Error merging BAM files: {e}")

# Run samtools sort
try:
    subprocess.run(["samtools", "sort", output_merged_bam, "-o", sorted_merged_bam], check=True)
    print(f"Sorted merged BAM file saved as {sorted_merged_bam}")
except subprocess.CalledProcessError as e:
    print(f"Error sorting merged BAM file: {e}")

# Define file names for PCR duplicates removal
input_MergeX = "MergeInsertionsS.bam"
output_MergeX = "MergeInsertionsD.bam"

# Check if the sorted BAM file exists before proceeding
if not os.path.exists(input_MergeX):
    raise FileNotFoundError(f"Sorted BAM file not found: {input_MergeX}")

# Run samtools rmdup (remove PCR duplicates)
try:
    subprocess.run(["samtools", "rmdup", "-sS", input_MergeX, output_MergeX], check=True)
    print(f"Removed duplicates and saved as {output_MergeX}")
except subprocess.CalledProcessError as e:
    print(f"Error removing duplicates: {e}")

# Run samtools index
try:
    subprocess.run(["samtools", "index", output_MergeX], check=True)
    print(f"Indexed BAM file: {output_MergeX}")
except subprocess.CalledProcessError as e:
    print(f"Error indexing BAM file: {e}")



#---------------------------------------------------------------------------------------------------------------


##filter for soft clipped reads, make .fai file from Genome
# Define file names
input_MergedY = "MergeInsertionsS.bam"
output_Clipped = "Clipped.bam"

# Check if input BAM file exists
if not os.path.exists(input_MergedY):
    raise FileNotFoundError(f"Input BAM file not found: {input_MergedY}")

# Run samtools view with awk filtering for soft-clipped reads
try:
    # Filter for soft-clipped reads
    awk_command = "awk '$6 ~ /S/{print}; $1 ~ /@/{print}'"
    samtools_view_cmd = f"samtools view -h {input_MergedY} | {awk_command} | samtools view -bS - > {output_Clipped}"
    
    subprocess.run(samtools_view_cmd, shell=True, check=True)
    print(f"Filtered soft-clipped reads and saved as {output_Clipped}")
except subprocess.CalledProcessError as e:
    print(f"Error filtering soft-clipped reads: {e}")

# Index the output BAM file
try:
    subprocess.run(["samtools", "index", output_Clipped], check=True)
    print(f"Indexed BAM file: {output_Clipped}")
except subprocess.CalledProcessError as e:
    print(f"Error indexing BAM file: {e}")

# Define source and destination paths for file copy
source = "/home/mga/Seq_data2/Clipped.bam"
destination = "/home/mga/Seq_data1/Clipped.bam"

# Ensure the destination directory exists
destination_dir = os.path.dirname(destination)
if not os.path.exists(destination_dir):
    os.makedirs(destination_dir)

# Copy the file
try:
    shutil.copy2(source, destination)
    print(f"Copied {source} to {destination}")
except Exception as e:
    print(f"Error copying the file: {e}")

# Define input and output file names
input_fasta = "Arab_Genome.fa"

# Run samtools faidx to create .fai index file
subprocess.run(["samtools", "faidx", input_fasta], check=True)

print(f"FASTA index file created: {input_fasta}.fai")



#--------------------------------------------------------------------------------------------------------------




#Update git repository
# Define the repository path
repo_path = "/home/mga/Seq_data1"

# Step 1: Check if there are any changes to commit
status = subprocess.run(["git", "status", "--porcelain"], cwd=repo_path, check=True, text=True, capture_output=True)

if status.stdout.strip():  # If there are any changes (not empty output)
    # Step 2: Add changes
    subprocess.run(["git", "add", "."], cwd=repo_path, check=True)
    
    # Step 3: Commit changes with a message
    subprocess.run(["git", "commit", "-m", "Adding Genome Files"], cwd=repo_path, check=True)
    
    # Step 4: Push changes to the remote repository
    subprocess.run(["git", "push", "origin", "main"], cwd=repo_path, check=True)
    
    print("Git repository updated successfully!")
else:
    print("No changes to commit.")
