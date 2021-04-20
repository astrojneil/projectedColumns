#Projected Column Densities
Highlights of analysis code to generate projected column densities with Trident for various UV backgrounds generated with Cloudy


###To generate new background tables for Trident to use: 
1. Download Cloudy
2. Put the folders 'Cooling_tools' and 'HM_shu' into the top level of the Cloudy folders
3. Inside HM_shu; edit hm_2012.par
  - edit the cloudy path to your cloudy executable
  - edit the output directory and prefix for where you want the files to go/what they should be called
  - edit the line "command table HM12 redshift 0.0 [factor=10000.00] to the correct parameters for your background
4. Inside HM_shu; run: 
  $ ./CIAOLoop hm_2012.par      (you can include -np 8 if you have multiple cores)
  
5. Inside your output directory; 
 - copy convert.py and cloudy_ascii_hdf5.py to here
 - run convery.py with the following
  $ python convert.py *name of .run file*  *name of final table*  H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V CR Mn Fe Co Ni Cu Zn
  
  
###To generate column densitiy projections:
1. Run makeColDens_db.py

Inside the file, you will need to change; 
  - the paths to the simulation chk files
  - the names of the database/output files
  - the path to the ionization table generated with Cloudy
  
