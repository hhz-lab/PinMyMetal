#!/usr/bin/python
import getopt
import os, sys
import subprocess

script_dir = os.path.dirname(os.path.abspath(__file__))

options, remainder = getopt.getopt(sys.argv[1:], 'i:o', ['intput='])
for opt, arg in options:
    if opt in ('-i', '--input'):
        inputid = arg

pdbid=inputid

# Define the database name
dbname = "metal" + str(pdbid)

data_dir = os.path.join(script_dir, f"{pdbid}_nb_result")

def execute_with_psql(dbname, sqlfile):
    """
    Executes an SQL script using the `psql` command-line tool.

    Parameters:
        dbname (str): The name of the database.
        sqlfile (str): The path to the SQL script.
    """
    try:
        # Construct the `psql` command
        command = [
            "psql",
            "-h", "postgresql_container",  # Database host (container name)
            "-U", "myuser",                # Database user
            "-d", dbname,                  # Target database
            "-p", "5432",                  # Port number
            "-f", sqlfile                  # Path to the SQL script
        ]

        # Execute the command
        env = os.environ.copy()
        env["PGPASSWORD"] = "123456"  # Pass the password via environment variable
        result = subprocess.run(command, env=env, check=True, text=True, capture_output=True)

        # Print success message
        print(f"Executed {sqlfile} successfully.")
        print(result.stdout)

    except subprocess.CalledProcessError as e:
        print(f"Error executing {sqlfile}: {e.stderr}")

# Execute the scripts
execute_with_psql(dbname, os.path.join(script_dir, "createNeighborhoodPartialTables.sql"))
execute_with_psql(dbname, os.path.join(data_dir, "copyNeighborhoodData.sql"))

