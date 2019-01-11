How to build:
    
To build the program, you need CMake and Eigen installed. The CMake script and make script are shown as below:
       
$ mkdir /path/to/build
$ cd /path/to/build
$ cmake /path/to/src
$ make
 
run script:
$ cd /path/to/build
$ cd bin
$ ./main "/path/to/matrix/L" "sparse/dense" "/path/to/vector/b"

You need to specify the right-hand side vector is whether sparse or dense for now.
