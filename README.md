# etching_featue_profile
1. Create your local repository
   mkdir ./etching_feature_profile
2. Set remote origin
   git remote add origin https://github.com/hawkingeric/etching_feature_profile.git
3. Pull this repository (master branch) to your local site 
   git pull origin master
4. Enter single_finite_trench/single_long_trench/array_set directory or create your own
   cd single_finite_trench  or
   cd single_long_trench    or
   cd array_set
5. Compile with g++  
   g++ ../\*.cpp  
   g++ -o etching_openMP.o -fopenmp ../\*.cpp
6. JSON file "json_test.json" contains the input setting. Please edit it appropriately for your own purpose.
   1. Specify your output file directory under tag "output_file_directory".
   2. Specify your IEADF_data directory in json file under tag "IEADF_data_directory".
7. Run the job
   ./etching_openMP.o -input_json json_test.json
