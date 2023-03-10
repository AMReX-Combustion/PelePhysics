template <typename T>
void print(T * array, int size, std::string name="") {
    if (name!="") {
        std::cout << name << " : "  << "\n";
    }
    for (int i=0; i<size; i++) {
        if (std::abs(array[i])>1e-300){
            std::cout << "\t" << i << " : " <<  array[i] << "\n";
        }
    }
}

template <typename T>
void print(T val,  std::string name="") {
    if (name!="") {
        std::cout << name << " : "  << "\n";
    }
    std::cout << "\t" <<  val << "\n";
}
