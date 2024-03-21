# Instructions to run the DBEM

Before you run the DBEM make sure you read this document, otherwise things may not function properly and your work will be delayed. As a overall suggestion, I recommend you use `git` and sync it with Compute Canada. That way you can modify everything on your local computer and then simply `push + pull` to Compute Canada.

## Settings file

The `settings.txt` file contains the main information to run the DBEM. Make sure all the info is correct:

Table 1: Glossary for settings file

| Vector    | Description                                                                                                                                                    | Example         |
|----------------|----------------------------------------|----------------|
| SppNo     | Number of species to run                                                                                                                                       | 40              |
| CCSc      | Climate change data                                                                                                                                            | C6GFDL26        |
| SSP       | SSP data for fishing effort. Only applicable if fishing effort is activated.                                                                                   | SSP126          |
| rsfile    | Name of file containing the list of species to run                                                                                                             | RunSppList      |
| rpath     | Folders where results will be created in your `scratch` folder (auto. by DBEM)                                                                                 | C6GFDL26F1MPA10 |
| tpath     | Folder containing taxon information                                                                                                                            | TaxonDataC0     |
| ifile     | Number of species list to run. Note this number will automatically increase every time you run the DBEM                                                        | 10              |
| FHS       | Fishing level in the high seas as a proportion of MSY                                                                                                          | 1.00            |
| FEEZ 1.00 | Fishing level in EEZs as a proportion of MSY                                                                                                                   | 0.50            |
| MPApath   | Name of the MPA scenario to run. Note: The name of the vector need to match the file name in the `Data` folder. If you don't have an MPA scenario use `mpa_no` | mpa_30          |

\
\

## Running file

You will need to prepare the *slurm* process to run the DBEM using the `run_dbem.sh` script. At the minimum, you need to modify the following variables before running the script:

Table 2 : Minimum information needed to run the DBEM. See [here](https://hpc.nmsu.edu/discovery/slurm/commands/) if you want to know what each vector means and how include more/less.

| Vector      | Description                                         | Example (standard)                                                    |
|-----------------|----------------------|---------------------------------|
| N           | Nodes needed                                        | 1                                                                     |
| mem-per-cpu | Memory needed for the DBEM                          | 700M                                                                  |
| t           | Computing time dd-hh:mm:ss                          | 01-12:00:00                                                           |
| mail-user   | Email account where you will get run notificationds | youremail[\@oceans.ubc.ca](mailto:mail-user=j.palacios@oceans.ubc.ca) |

**Note** that you do not need to modify anything after line 11. However, if you changed the name or location of your source script, then you will need to update that on line 20

``` bash
# This is the path that runs the DBEM 
../../dbem_scripts/DBEM_v2_y
```

## Compile Code

Before you actually run the DBEM source code **for the first time** , you need to compile the Fortran code (`DBEM_v2_y.f90` or `DBEM_v2_m.f90`). You will also need to compile the code every time you make a modification to the source code. For that, you can simply run the `compile_dbem.sh` script in the *terminal* of Cedar at *Compute Canada* as follows:

``` bash
# First we give permitions to execute the script
chmod +x compile_dbem.sh

# Now we execute the script
./compile_dbem.sh
```

You should get *no return message* but if you write `ls` in the terminal you should see new files (e.g., *dbem.mod*). Every time you modify the DBEM original.

## Run the DBEM

Once you complete all of those steps you are ready to run the DBEM in the *terminal* of Cedar at *Compute Canada* as follows:

``` bash

emacs run_dbem.sh
```

