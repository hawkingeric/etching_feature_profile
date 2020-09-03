# etching_featue_profile
1. Download the binary execution file in bin/etching (without openMP) or etching_openMP.o (with openMP).
   Or you could download all the .cpp .h .cpp files and compile by yourself.
   g++ *.cpp
   g++ -o etching_openMP.o -fopenmp *.cpp
2. Download json_test.json.
3. Specify your output file directory under tag "output_file_directory".
4. Download IEADF_data and specify your IEADF_data directory in json file under tag "IEADF_data_directory".
