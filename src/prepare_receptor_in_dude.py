import os
import subprocess

def process_pdb_files_and_run_command(root_path):
    """
    Process all receptor.pdb files in the given root directory.
    Deletes lines containing specific substrings and replaces the 13th character with a space if it's not a space.
    After processing each file, it changes the directory to the file's directory and runs a command.
    If the command returns an error, it prints the error message.
    """
    for subdir, dirs, files in os.walk(root_path):
        for file in files:
            if file == "receptor.pdb":
                file_path = os.path.join(subdir, file)

                # Process the file
                with open(file_path, 'r') as f:
                    lines = f.readlines()

                with open(file_path, 'w') as f:
                    for line in lines:
                        # Check for unwanted substrings and skip the line if any are found
                        if any(substr in line for substr in ["ZN","MG"]):
                            continue

                        # Replace 13th character with a space if it's not a space
                        if line.startswith("ATOM")  and line[12] != ' ': 
                            line = line[:12] + ' ' + line[13:]

                        f.write(line)

                # Change directory and run the command
                os.chdir(subdir)
                command = "prepare_receptor -r receptor.pdb -A hydrogens  "
                try:
                    subprocess.run(command, shell=True, check=True)
                except subprocess.CalledProcessError as e:
                    print(f"An error occurred while running command in {subdir}: {e}")

# Define the root path where the directories are located
root_path = '/mnt/e/SDT/input/all'

# Process the pdb files and run the command for each
process_pdb_files_and_run_command(root_path)

print("Processing and command execution complete.")

