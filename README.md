# RDF-ntuple-Cpp
RDataFrame C++ version of ntuple software

C++ version of ntuple software

To run, have all .h and .cxx files in the same directory, then run:



 - `rootcling -v4 -f mydict.cxx  -rmf libmydict.rootmap -rml libmydict.so  LinkDef.h`



 - `g++ -shared -o libmydict.so mydict.cxx `\``root-config --cflags --libs`\`` -fPIC`



After the dictionary is created for vector vector, then compile

 - `g++ -g -o RDF_analyse_h_Za RDF_analyse_h_Za.cxx MakeRDF.cxx RDFObjects.cxx RDFevent.cxx -O2 $(root-config --libs --cflags) -std=c++14`


And run


`./RDF_analyse_h_Za`

