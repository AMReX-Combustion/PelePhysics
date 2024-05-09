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
