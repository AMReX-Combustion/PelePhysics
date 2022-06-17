template <typename T>
void print(T * array, int size, std::string name="") {
    if (name!="") {
        std::cout << name << " : "  << "\n";
    }
    for (int i=0; i<size; i++) {
        std::cout << "\t" << i << " : " <<  array[i] << "\n";
    }
}

template <typename T>
void print(T val,  std::string name="") {
    if (name!="") {
        std::cout << name << " : "  << "\n";
    }
    std::cout << "\t" <<  val << "\n";
}

template <typename T>
void print(std::string filename, T * array, int size, std::string name="") {
    std::ofstream file;
    file.open(filename,std::ios_base::app);
    

    if (name!="") {
        file << name << " : "  << "\n";
    }
    for (int i=0; i<size; i++) {
        file << "\t" << i << " : " <<  array[i] << "\n";
    }
    file.close();
}

template <typename T>
void print(std::string filename, T val,  std::string name="") {
    std::ofstream file;
    file.open(filename, std::ios_base::app);
 
    if (name!="") {
        file << name << " : "  << "\n";
    }
    file << "\t" <<  val << "\n";
    file.close();
}

template <typename T>
void print2py(std::string filename, T * array, int size, std::string name="") {
    std::ofstream file;
    file.open(filename,std::ios_base::app);

    file << "import numpy as np"  <<  "\n";
    file << name << " = np.zeros("  << size <<  ")\n";

    for (int i=0; i<size; i++) {
        file << name << "[" << i << "]=" <<  array[i] << "\n";
    }
    file.close();
}
