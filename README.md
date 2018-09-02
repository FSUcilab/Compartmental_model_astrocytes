# Compartmental_model_astrocytes
A simple compartmental model for astrocytes to explore interactions between major processes and soma

# Compartmental_model_astrocytes
A simple compartmental model for astrocytes to explore interactions between major processes and soma


C code to explore the interaction between several major processes in an astrocyte and a soma. The code allows duplication of results found in the paper 
"A Compartmental Model to Investigate Local and Global Calcium Dynamics in Astrocytes" submitted to Frontiers. 

To run the code: 

      g++ -o astro -c astro_2017.c normal.c
      
      
To execute the code: 

   ./astro time nsteps der sig dc corr addOn out_file seed
   
where: 

   time:     duration of simulation in non-dimensional time unites
   nsteps:   total number of time steps
   der:      diffusivity in the ER
   sig:      sigma: standard deviation of neuronal input
   dc:       diffusivity in the cytosol
   corr:     pairwise correlation coefficient
   addOn:    
   out_file: name of output file that contains time signals
   see:      initial seed for random input signals
   
   ===================================
   
   Python 2.7 is required.
   
   Alternatively, run the python script: 
   
      python create_multiple_runs.py
      
   which will take care of compileing the code and running it. To run multiple runs, set the "runs" variable
   in create_multiple_runs.py . 
   
   Update the following variables in input.py to set the proper paths: 
   
        base = os.path.join(repo_dir, "comp_astro_evan", "refactored")
        libs = os.path.join(repo_dir, base, "libs")
        
  =============================================================================
  
  Once create_multiple_runs.py has completed, directories run_data000/, run_data001/, ... will have been created, 
  which will contain output files, and csv files with spiking data. 
  
