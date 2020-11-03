# etching_featue_profile
1. Compile with g++
   g++ *.cpp
   g++ -o etching_openMP.o -fopenmp *.cpp
2. JSON file "json_test.json" contains the input setting. Please edit it appropriately for your own purpose.
3. Specify your output file directory under tag "output_file_directory".
4. Directory "IEADF_data" contains the ion-energy-angular distribution function as needed.
5. Specify your IEADF_data directory in json file under tag "IEADF_data_directory".
6. Run the job
   ./{executable file} -input_json json_test.json
