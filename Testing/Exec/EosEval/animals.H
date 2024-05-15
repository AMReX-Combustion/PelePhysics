#ifndef ANIMALS_H
#define ANIMALS_H

#include "Factory.H"

struct lab;
struct doberman;
struct pitbull;


namespace pele::physics {

class animal
{
  int kingdom = 0;
};

  class data_cont {
  public:
    virtual void dcont_type() {amrex::Print() << "dcont type 0" << std::endl;}
  };
  class data_cont1 : data_cont {
  public:
    void  dcont_type() override { amrex::Print() << "dcont type 1" << std::endl;}
    };
  
  class data_cont2 : data_cont {
  public:
    void  dcont_type() override { amrex::Print() << "dcont type 2" << std::endl;}
    };
  
  template <typename CatsType>
  class CATS : public Factory<CATS<CatsType>>
  {
  public: 
    static std::string base_identifier(){ return "base_cat_arbitrary";}

    virtual void print_hello() {amrex::Print() << "HELLO from base arbitrary" <<std::endl;}
    
    virtual std::string name(); // { return "I am Fred the base cat";}
    
    CatsType my_cont;
    
    virtual void cont_name() {my_cont.dcont_type();}
    
  };

  /*
  template <>
  class CATS<data_cont> : public Factory<CATS<data_cont>>
  {
  public: 
    static std::string base_identifier(){ return "base_cat_data_cont";}

    virtual std::string name(); // { return "I am Fred the base cat";}
  };
  */
  
  class Smagorinsky : public CATS<data_cont>::Register<Smagorinsky> {
  public:
    static std::string identifier() { return "Smagorinsky"; }
    std::string name() override; // { return "I am Fred the Smag cat";}
  };
  
  class Lion : public CATS<data_cont>::Register<Lion> {
  public:
    static std::string identifier() { return "Lion"; }
    std::string name() override;// { return "I am Fred the Lion";}
  };


  template <>
  class CATS<data_cont2> : public CATS<data_cont>::Register<CATS<data_cont2>> {
  public:
    static std::string identifier() { return "Cats2"; }
    std::string name() override;// { return "I am Fred the Lion";}
    data_cont2 my_cont;
  };
  
  class LIZARDS : public Factory<LIZARDS>
  {
  public: 
    static std::string base_identifier(){ return "base_lizard";}
    virtual std::string name(); // { return "I am Fred the generic lizard";}
  };
  
  class geico : public LIZARDS::Register<geico> {
  public:
    static std::string identifier() { return "geico"; }
    std::string name() override; // { return "I am Fred the Geico";}
  };
  
  class croc : public LIZARDS::Register<croc> {
  public:
    static std::string identifier() { return "croc"; }
    std::string name() override; // { return "I am Fred the Croc";}
  };

  std::string LIZARDS::name() { return "I am Fred the generic";}
  std::string croc::name() { return "I am Fred the croc";}
  std::string geico::name() { return "I am Fred the Geico";};

  class dog_type { public: void type() {amrex::Print() << "type is GENERIC" << std::endl;};};
  class dog_doberman : public dog_type {public: void type() {amrex::Print() << "type is DOBERMAN" << std::endl;};};
  class dog_golden : public dog_type { public: void type() {amrex::Print() << "type is GOLDEN" << std::endl;};};
  
  template <typename DogType>
  class dog {
  public:
    virtual void type() { std::cout << "IN FUNC "; my_type.type();};
    DogType my_type;
  };

  /*
  template <typename DogType, typename BaseType>
  class dogplus : public dog<DogType> {
  public:
    virtual void type() { std::cout << "IN FUNC "; my_type.type();};
    DogType my_type;
    }; */

  class data_base {
  public:
    virtual void print() {std::cout  << "I am data_base" <<std::endl;}
    int x = 0;
  };

    class data_der1 : public data_base{
    public:
      void print() override {std::cout  << " I am data_der1" <<std::endl; x=7;}
      int y = 2;
    };
  
    class data_der2 : public data_base{
    public:
      void print() override {std::cout  << "I am data_der2" <<std::endl; x=23;}
      int z = 3;
    };
  
  template <typename BaseDataType>
  class base_class : public Factory<base_class<BaseDataType>>
  {
  public:
    static std::string base_identifier(){ return "base_class";}
    std::string name() {return "no name";};
    virtual BaseDataType* get_my_data() = 0;
    
  };

  template <typename ClassType, typename BaseDataType = ClassType>
  class derived_class : public  base_class<BaseDataType> //::template Register<derived_class<ClassType, BaseDataType>> 
    {
    public:
      static_assert(std::is_base_of_v<BaseDataType,ClassType> == true);
      std::string name()  {return "derived names";};
      // static std::string identifier() {return "derived_class";};
      virtual BaseDataType* get_my_data() {return &myData;};
    private:
      ClassType myData;
    };

  class fart : public base_class<data_base>::Register<fart, derived_class<data_der1, data_base>> {
    //class fart : public derived_class<data_der1, data_base> {
  public:
    std::string name();
    static std::string identifier() {return "fart_class";};
    //data_base* get_my_data() override {return &myData;};
  };
  
  class fart2 : public  base_class<data_base>::Register<fart2, derived_class<data_der2, data_base>> {
  //class fart2 : public  derived_class<data_der2, data_base> {
  public:
    std::string name();
    static std::string identifier() {return "fart2_class";};
    //data_base* get_my_data() override {return &myData;};
  };
  
  //std::string fart::name() {return "my name is fart class";};
  //std::string fart::identifier() {return "fart_class";};
  std::string fart::name () {return "my name is fart class";};
  //data_base* fart::get_my_data() {return nullptr;};
  //std::string fart2::identifier() {return "fart2_class";};
  std::string fart2::name () {return "my name is fart2 class";};
  //data_base* fart2::get_my_data() {return nullptr;};

  // fart a;
  // fart2 b; 

  /*
  class other_class : public base_class<data_base>::Register<other_class>
  {
  public:
    static std::string identifier();
    std::string name() override;
  };
  
  std::string other_class::name() {return "my name is other class";};
  std::string other_class::identifier() {return "other_class";};
  */

  /*

  template<>
  std::string derived_class<data_cont>::name() {return "my name is derived data_cont";}

  template<>
  std::string derived_class<data_cont>::identifier() {return "derived_class";}
  */
  
  class DOGS : public dog<dog_type>, public Factory<DOGS>
  {
  public: 
    static std::string base_identifier(){ return "generic_dog";}
    virtual std::string name();
  };
  
  class doberman : public dog<dog_doberman>, public DOGS::Register<doberman> {
  public:
    static std::string identifier() { return "doberman"; }
    std::string name() override;
  };
  
  class golden : public dog<dog_golden>, public DOGS::Register<golden> {
  public:
    static std::string identifier() { return "golden"; }
    std::string name() override; 
  };

  std::string DOGS::name() { return "I am Fred the generic";}
  std::string doberman::name() { return "I am Fred the doberdog";}
  std::string golden::name() { return "I am Fred the golden retriever";};

  template<>
  std::string CATS<data_cont>::name() { return "I am Fred the base cat";}
  
  std::string Smagorinsky::name()  { return "I am Fred the Smag cat";}
  std::string Lion::name() { return "I am Fred the Lion";}
  std::string CATS<data_cont2>::name() { return "I am Fred the 2nd cat";}
}

#endif
