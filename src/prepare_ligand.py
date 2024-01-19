import os 
import subprocess

def process_pdb_files_and_run_command(root_path):
   
    for subdir, dirs, files in os.walk(root_path):
        for file in files:
            if file == "receptor.pdb":
                file_path = os.path.join(subdir, file)

                

                # Change directory and run the command
                os.chdir(subdir)
                command = "prepare_receptor -r receptor.pdb -A hydrogens  "
                try:
                    subprocess.run(command, shell=True, check=True)
                except subprocess.CalledProcessError as e:
                    print(f"An error occurred while running command in {subdir}: {e}")

